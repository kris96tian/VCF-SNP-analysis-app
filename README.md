# VCF/SNP Analyzer App

Link to the Web-App: https://vcfanalysis.pythonanywhere.com/

A Flask web application for analyzing VCF (Variant Call Format) files containing genetic variants such as SNPs (Single Nucleotide Polymorphisms).

<img src="blob:chrome-untrusted://media-app/c61587a6-fbd4-473d-8e6f-38bfd86 a6080"/>![image](https://github.com/user-attachments/assets/2593f7c5-b05b-4be5-b588-ce9a6e770779)


## Features

- Upload and validate VCF files (supports .vcf, .vcf.gz, and .bcf formats)
- Automatic fixing of common VCF format issues
- Comprehensive variant analysis including:
  - Variant type distribution
  - Chromosome distribution
  - Quality score distribution
  - Transition/Transversion ratio
  - Population genetics metrics
- Interactive visualizations using Plotly
- Asynchronous processing of large files with job status tracking

## Installation

### Prerequisites

- Python 3.8+
- pip

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/snp-analyzer-app.git
   cd snp-analyzer-app
   ```

2. Create and activate a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Start the application:
   ```bash
   python app.py
   ```

2. Open your web browser and navigate to `http://127.0.0.1:5000/`

3. Upload a VCF file:
   - Click on "Choose File" and select a VCF file
   - Click "Analyze File" to process the file and view results
   - Alternatively, click "Validate Only" to check if the file is valid without running analysis

4. View analysis results:
   - Chromosome distribution visualization
   - Quality score distribution
   - Variant type statistics
   - Ts/Tv ratio and other metrics

## Technical Details

- **Backend**: Flask (Python)
- **Data Processing**: cyvcf2, numpy
- **Frontend**: HTML, JavaScript, Tailwind CSS
- **Visualization**: Plotly
- **File Management**: Only temporary files used during processing (no persistent storage)

## Development

To contribute to this project:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin feature-name`
5. Submit a pull request
