import os
import click
import zipfile
import xml.etree.ElementTree as ET

# it assumes that the template xml file is in the same directory as the script
script_dir = os.path.dirname(__file__)
template_file = os.path.join(script_dir, "sampleOnlyUpload.xml")  # The original XML template

def create_xml_fun(template, output_xml, vcf_files):
    """
    Uses a template XML file to create a XML file with sample information and variant filenames,
    and saves the XML to a specified output file.
    Args:
        template (str): Path to the template XML template file.
        output_xml (str): Path to save the modified XML file.
        vcf_files (list of str): List of VCF file paths to extract sample names and filenames.
    Functionality:
        - Parses the input XML template file.
        - Updates the <Name> tag in the XML with the sample name derived from the first VCF file.
        - Adds or updates the <SubjectId> tag with the same sample name.
        - Clears existing <VariantsFilenames> entries and adds new entries for each VCF file.
        - Ensures proper namespace handling and formatting in the output XML.
    Notes:
        - The sample name is derived from the first VCF file in the `vcf_files` list.
        - The function ensures that the XML output adheres to the namespace defined in the template.
    Raises:
        FileNotFoundError: If the template file does not exist.
        ET.ParseError: If the XML template file is not well-formed.
    """

    tree = ET.parse(template)
    root = tree.getroot()

    # Define the namespace
    ns = {"qci": "http://qci.qiagen.com/xsd/interpret"}

    # Register the namespace to avoid ns0 prefix in the output_xml
    ET.register_namespace('', ns["qci"])

    # Update the <Name> tag
    name_element = root.find(".//qci:Sample/qci:Name", ns)
    if name_element is not None:
        for vcf in vcf_files:
            sample = os.path.splitext(os.path.splitext(os.path.basename(vcf))[0])[0] + "_TEST" #TODO: remove test
            break
        name_element.text = sample

        # Add SubjectId after Name and before VariantsFilenames
        subject_id_element = root.find(".//qci:Sample/qci:SubjectId", ns)
        if subject_id_element is None:
            name_element = root.find(".//qci:Sample/qci:Name", ns)
            if name_element is not None:
                # Insert SubjectId after Name
                sample_element = root.find(".//qci:Sample", ns)
                subject_id_element = ET.Element("SubjectId")
                subject_id_element.text = sample
                subject_id_element.tail = "\n"
                name_index = list(sample_element).index(name_element)
                sample_element.insert(name_index + 1, subject_id_element)
        else:
            subject_id_element.text = sample
            subject_id_element.tail = "\n"
        
    # Find the <VariantsFilenames> element and clear existing filenames
    filenames_element = root.find(".//qci:Sample/qci:VariantsFilenames", ns)
    if filenames_element is not None:
        # Remove existing children (if any)
        for child in list(filenames_element):
            filenames_element.remove(child)

        # Add new filenames
        for filename in vcf_files:
            new_file_element = ET.SubElement(filenames_element, "VariantsFilename")
            new_file_element.text = os.path.basename(filename)
            # Add a newline before each new element
            new_file_element.tail = "\n"

    # Save the modified XML to a new file
    tree.write(output_xml, encoding="UTF-8", xml_declaration=True)

def create_zip_fun(xml_file, vcf_files, zip_file):
    """
    Creates a ZIP archive containing an XML file and VCF files.
    Args:
        xml_file (str): The path to the XML file to be added to the ZIP archive.
        vcf_files (list of str): A list of paths to VCF files to be added to the ZIP archive.
        zip_file (str): The path where the resulting ZIP archive will be created.
    Returns:
        None
    """
    
    with zipfile.ZipFile(zip_file, 'w') as zipf:
        # Add the XML file to the zip
        zipf.write(xml_file, os.path.basename(xml_file))
        
        # Add each VCF file to the zip
        for vcf in vcf_files:
            zipf.write(vcf, os.path.basename(vcf))

def send_qci_api_fun(zip_file, api_key_file):
    """
    Sends a ZIP file to the QCI API using a provided API key.
    This function reads an API key from a specified file and uses it to 
    authenticate a POST request to the QCI API. The request uploads a ZIP 
    file to the API endpoint.
    Args:
        zip_file (str): The path to the ZIP file to be uploaded.
        api_key_file (str): The path to the file containing the API key.
    Raises:
        Exception: If an error occurs during the execution of the API command.
    Notes:
        - The API key file should contain the API key as a single line of text.
        - The function uses the `curl` command-line tool to send the request.
        - Ensure that the `curl` tool is installed and available in the system's PATH.
    """
    
    with open(api_key_file, 'r') as key_file:
        api_key = key_file.read().strip()
        
    try:
        api_command = f'curl -X POST "https://api.qiagenbioinformatics.eu/v2/sample" -H "accept: application/json" -H "Authorization: {api_key}" -H "Content-Type: multipart/form-data" -F "file=@{zip_file};type=application/zip"'
        os.popen(api_command).read()

    except Exception as e:
        print(f"An error occurred: {e}")
    
def complete_qci_submission_fun(template, output_xml, vcf_files, zip_file, api_key_file):
    """
    Core logic for creating XML, zipping files, and submitting to QCI.
    """
    try:
        create_xml_fun(template, output_xml, vcf_files)
    except Exception as e:
        raise RuntimeError(f"Error during XML creation: {e}")

    try:
        create_zip_fun(output_xml, vcf_files, zip_file)
    except Exception as e:
        raise RuntimeError(f"Error during ZIP creation: {e}")

    try:
        send_qci_api_fun(zip_file, api_key_file)
    except Exception as e:
        raise RuntimeError(f"Error during API submission: {e}")

@click.group(help="Create the QCI XML file for a sample. Submit the API command to QCI to upload the zip with xml and vcf files.")
def cli():
    pass

@cli.command(help="Create a QCI XML file for a sample.")
@click.option('--template', type=click.Path(exists=True), default=template_file, required=True, help='Path to the XML template file.', show_default=True)
@click.option('--output_xml', type=click.Path(), required=True, help='Path to the output XML file.', show_default=True)
@click.option('--vcf_files', type=click.Path(exists=True), multiple=True, required=True, help='List of variant filenames to be added to the XML. Example: --vcf_files file1.vcf --vcf_files file2.vcf', show_default=True)
def create_xml(template, output_xml, vcf_files):
    create_xml_fun(template, output_xml, vcf_files)

@cli.command(help="Create a zip file with an XML file and multiple VCF files.")
@click.option('--xml_file', type=click.Path(exists=True), required=True, help='Path to the XML file to be included in the zip.')
@click.option('--vcf_files', type=click.Path(exists=True), multiple=True, required=True, help='List of VCF filenames to be included in the zip. Example: --vcf_files file1.vcf --vcf_files file2.vcf')
@click.option('--zip_file',type=click.Path(),  required=True, help='Path to the output zip file.')
def create_zip(xml_file, vcf_files, zip_file):
    create_zip_fun(xml_file, vcf_files, zip_file)

@cli.command(help="Send the QCI API call to upload the zip file.")
@click.option('--zip_file',type=click.Path(), required=True, help='Path to the output zip file.', show_default=True)
@click.option('--api_key_file', type=click.Path(exists=True), required=True, help='File with QCI API key.', default = "/clinical/exec/wgs_somatic/dependencies/credentials/qci_api_key.txt", show_default=True)
def submit_sample_qci(zip_file, api_key_file):
    send_qci_api_fun(zip_file, api_key_file)

@cli.command(help="All-in-one command to submit a sample (upload vcf file) to QCI. Create XML, zip it with vcf file, and submit to QCI with API.")
@click.option('--template', type=click.Path(exists=True),default=template_file, required=True, help='Path to the XML template file.', show_default=True)
@click.option('--output_xml', type=click.Path(), required=True, help='Path to the output XML file.', show_default=True)
@click.option('--vcf_files', type=click.Path(exists=True), multiple=True, required=True, help='List of variant filenames to be added to the XML. Example: --vcf_files file1.vcf --vcf_files file2.vcf', show_default=True)
@click.option('--zip_file',type=click.Path(), required=True, help='Path to the output zip file.', show_default=True)
@click.option('--api_key_file', type=click.Path(exists=True), required=True, help='File with QCI API key.', default = "/clinical/exec/wgs_somatic/dependencies/credentials/qci_api_key.txt", show_default=True)
def complete_qci_submission(template, output_xml, vcf_files, zip_file, api_key_file):
    try:
        complete_qci_submission_fun(template, output_xml, vcf_files, zip_file, api_key_file)
    except RuntimeError as e:
        print(e)
        
if __name__ == "__main__":
    cli()
    