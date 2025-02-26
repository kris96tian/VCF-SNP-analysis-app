from flask import Flask, request, jsonify, render_template, send_from_directory, url_for
import os
import subprocess
import tempfile
from werkzeug.utils import secure_filename
import io
import base64
import numpy as np
from collections import defaultdict
import cyvcf2
import plotly.graph_objects as go
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
import logging
from pathlib import Path
import uuid

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config.from_mapping(
    UPLOAD_FOLDER='uploads',
    OUTPUT_FOLDER='outputs',
    ALLOWED_EXTENSIONS={'vcf', 'vcf.gz', 'bcf'},
    MAX_CONTENT_LENGTH=500 * 1024 * 1024,
    PROCESS_EXECUTOR=ThreadPoolExecutor(max_workers=4)
)


Path(app.config['UPLOAD_FOLDER']).mkdir(exist_ok=True)
Path(app.config['OUTPUT_FOLDER']).mkdir(exist_ok=True)

# Job tracking dictionary
job_tracker = {}

class VCFPreprocessor:

    @staticmethod
    def validate_and_fix(filepath):
        fixed_filepath = str(filepath) + ".fixed"
        
        try:
            try:
                reader = cyvcf2.VCF(str(filepath))
                reader.close()
                return str(filepath), None  
            except Exception as e:
                logger.info(f"VCF validation failed, attempting to fix: {str(e)}")
                
            with open(filepath, 'r') as infile, open(fixed_filepath, 'w') as outfile:
                header_fixed = False
                for line in infile:
                    if line.startswith('#'):
                        if line.startswith('#CHROM'):
                            parts = line.strip().split()
                            fixed_line = '\t'.join(parts) + '\n'
                            outfile.write(fixed_line)
                            header_fixed = True
                        else:
                            outfile.write(line)
                    else:
                        if ' ' in line and '\t' not in line:
                            parts = line.strip().split()
                            fixed_line = '\t'.join(parts) + '\n'
                            outfile.write(fixed_line)
                        else:
                            outfile.write(line)
            
            try:
                reader = cyvcf2.VCF(fixed_filepath)
                reader.close()
                return fixed_filepath, None
            except Exception as e:
                return str(filepath), f"VCF validation failed after fixing attempt: {str(e)}"
                
        except Exception as e:
            return str(filepath), f"Error during VCF preprocessing: {str(e)}"

def allowed_file(filename):
    """Check if file extension is allowed"""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS'] or \
           filename.endswith('.vcf.gz')



class VCFAnalyzer:
    
    def __init__(self, filepath):
        self.filepath = filepath
        self.reader = cyvcf2.VCF(filepath)
        self.samples = self.reader.samples
        
    def comprehensive_analysis(self):
        """Perform comprehensive variant analysis with population genetics metrics"""
        stats = {
            'chromosomes': defaultdict(int),
            'qualities': [],
            'depths': [],
            'allele_freqs': [],
            'types': defaultdict(int),
            'ts_tv_ratio': {'transition': 0, 'transversion': 0},
            'hwe': [],
            'inbreeding_coeff': [],
            'consequences': defaultdict(int)
        }
        
        for variant in self.reader:
            # Chromosome distribution
            stats['chromosomes'][variant.CHROM] += 1
            
            # Variant type
            stats['types'][variant.var_type] += 1
            
            # Quality metrics
            if variant.QUAL:
                stats['qualities'].append(variant.QUAL)
            
            # Depth metrics
            if variant.INFO.get('DP'):
                stats['depths'].append(variant.INFO['DP'])
                
            # Allele frequency
            if variant.aaf:
                stats['allele_freqs'].append(variant.aaf)
                
            # Ts/Tv ratio
            if variant.is_snp:
                ref, alt = variant.REF, variant.ALT[0]
                if self._is_transition(ref, alt):
                    stats['ts_tv_ratio']['transition'] += 1
                else:
                    stats['ts_tv_ratio']['transversion'] += 1
                    
            # Population genetics metrics
            if 'HWE' in variant.INFO:
                stats['hwe'].append(variant.INFO['HWE'])
            if 'InbreedingCoeff' in variant.INFO:
                stats['inbreeding_coeff'].append(variant.INFO['InbreedingCoeff'])
                
            # Variant consequences
            if 'CSQ' in variant.INFO:
                for consequence in variant.INFO['CSQ'].split(','):
                    stats['consequences'][consequence.split('|')[1]] += 1
                    
        return self._process_stats(stats)
    
    def _is_transition(self, ref, alt):
        return (ref in 'AG' and alt in 'AG') or (ref in 'CT' and alt in 'CT')
    
    
    def _process_stats(self, stats): # derived statistics from raw data
        stats['ts_tv_ratio'] = (stats['ts_tv_ratio']['transition'] / 
                               stats['ts_tv_ratio']['transversion']) if stats['ts_tv_ratio']['transversion'] else 0
        stats['quality_stats'] = self._calculate_summary_stats(stats['qualities'])
        stats['depth_stats'] = self._calculate_summary_stats(stats['depths'])
        stats['af_stats'] = self._calculate_summary_stats(stats['allele_freqs'])
        stats['hwe_stats'] = self._calculate_summary_stats(stats['hwe'])
        stats['inbreeding_stats'] = self._calculate_summary_stats(stats['inbreeding_coeff'])
        return stats
    
    @staticmethod
    def _calculate_summary_stats(data):
        if not data:
            return {}
        arr = np.array(data)
        return {
            'mean': float(np.mean(arr)),
            'median': float(np.median(arr)),
            'std': float(np.std(arr)),
            'min': float(np.min(arr)),
            'max': float(np.max(arr)),
            'q25': float(np.quantile(arr, 0.25)),
            'q75': float(np.quantile(arr, 0.75))
        }

@app.route('/')
def home():
    return render_template('index.html')



@app.route('/api/upload', methods=['POST'])
def handle_upload():
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    file = request.files['file']
    if not file or file.filename == '':
        return jsonify({'error': 'Empty filename'}), 400
    
    if not allowed_file(file.filename):
        return jsonify({'error': 'Invalid file type'}), 400
    
    try:
        filename = secure_filename(file.filename)
        suffix = os.path.splitext(filename)[1]
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp_file:
            file.save(tmp_file.name)
            tmp_path = tmp_file.name
        job_id = str(uuid.uuid4())
        future = app.config['PROCESS_EXECUTOR'].submit(process_upload, tmp_path)
        job_tracker[job_id] = future
        
        return jsonify({
            'status': 'processing',
            'job_id': job_id,
            'results_url': url_for('show_results', filename=filename)
        }), 202
    except Exception as e:
        logger.error(f"Upload failed: {str(e)}")
        return jsonify({'error': 'File processing failed'}), 500



def process_upload(filepath):
    try:
        processed_filepath, validation_error = VCFPreprocessor.validate_and_fix(filepath)
        if validation_error:
            logger.error(f"VCF validation error: {validation_error}")
            return {"error": validation_error}
        analyzer = VCFAnalyzer(processed_filepath)
        result = analyzer.comprehensive_analysis()
        return result
    except Exception as e:
        logger.error(f"Analysis error: {str(e)}")
        return {"error": str(e)}
    finally:
        try:
            if os.path.exists(filepath):
                os.remove(filepath)
        except Exception as cleanup_error:
            logger.error(f"Error cleaning up uploaded file: {cleanup_error}")
        fixed_file = f"{filepath}.fixed"
        if os.path.exists(fixed_file):
            try:
                os.remove(fixed_file)
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up fixed file: {cleanup_error}")



@app.route('/results/<filename>')
def show_results(filename):
    job_id = request.args.get('job_id')
    return render_template('results.html', filename=filename, job_id=job_id)



@app.route('/api/results/<filename>')
def get_results(filename):
    job_id = request.args.get('job_id')
    if not job_id:
        return jsonify({'error': 'job_id required'}), 400
    if job_id in app.config.get('analysis_results', {}):
        results = app.config['analysis_results'][job_id]
        return jsonify({
            'status': 'complete',
            'results': results,
            'visualizations': generate_visualizations(results)
        })
    else:
        return jsonify({'status': 'processing'})



@app.route('/api/job/<job_id>')
def check_job_status(job_id):
    if job_id not in job_tracker:
        return jsonify({'error': 'Job not found'}), 404
    
    future = job_tracker[job_id]
    if future.done():
        try:
            result = future.result()
            if isinstance(result, dict) and 'error' in result:
                return jsonify({'status': 'failed', 'error': result['error']}), 500
            
            app.config.setdefault('analysis_results', {})
            app.config['analysis_results'][job_id] = result
            
            return jsonify({
                'status': 'complete',
                'results': result,
                'job_id': job_id
            })
        except Exception as e:
            return jsonify({'status': 'failed', 'error': str(e)}), 500
    else:
        return jsonify({'status': 'processing'})




def generate_visualizations(results):
    visuals = {}
    
    try:
        # Chromosome distribution
        if results.get('chromosomes'):
            chrom_fig = go.Figure()
            chrom_fig.add_trace(go.Bar(
                x=list(results['chromosomes'].keys()),
                y=list(results['chromosomes'].values()),
                marker_color='#6366f1'
            ))
            chrom_fig.update_layout(
                title='Variants by Chromosome',
                xaxis_title='Chromosome',
                yaxis_title='Count',
                template='plotly_dark'
            )
            visuals['chromosomes'] = chrom_fig.to_json()
        
        # quality distrib.
        if results.get('qualities'):
            qual_fig = go.Figure()
            qual_fig.add_trace(go.Histogram(
                x=results['qualities'],
                nbinsx=50,
                marker_color='#10b981'
            ))
            qual_fig.update_layout(
                title='Quality Score Distribution',
                xaxis_title='Quality',
                yaxis_title='Count',
                template='plotly_dark'
            )
            visuals['quality'] = qual_fig.to_json()
        
        return visuals
    except Exception as e:
        logger.error(f"Visualization generation error: {str(e)}")
        return {} 



@app.route('/api/validate', methods=['POST'])
def validate_vcf():
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    file = request.files['file']
    if not file or file.filename == '':
        return jsonify({'error': 'Empty filename'}), 400
    
    if not allowed_file(file.filename):
        return jsonify({'error': 'Invalid file type'}), 400
    
    try:
        filename = secure_filename(file.filename)
        suffix = os.path.splitext(filename)[1]
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp_file:
            file.save(tmp_file.name)
            tmp_path = tmp_file.name
        _, validation_error = VCFPreprocessor.validate_and_fix(tmp_path)
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception as cleanup_error:
            logger.error(f"Error cleaning up temporary file: {cleanup_error}")
        if validation_error:
            return jsonify({
                'valid': False,
                'error': validation_error
            })
        else:
            return jsonify({
                'valid': True,
                'message': 'VCF file is valid or was fixed successfully'
            })
    except Exception as e:
        logger.error(f"Validation failed: {str(e)}")
        return jsonify({'error': f'Validation failed: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(debug=True)