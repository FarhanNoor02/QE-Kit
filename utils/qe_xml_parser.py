# utils/qe_xml_parser.py
import xml.etree.ElementTree as ET
import os

def check_calc_status(prefix, outdir="./outdir/"):
    """
    Scans the outdir for prefix.xml and returns the last calculation type and status.
    In modern QE, the XML is usually at: {outdir}/{prefix}.save/data-file-schema.xml
    Or simply {prefix}.xml in older versions.
    """
    # Standard QE path for the XML data file
    xml_path = os.path.join(outdir, f"{prefix}.save", "data-file-schema.xml")
    
    # Fallback for local XML if not in the save folder
    if not os.path.exists(xml_path):
        xml_path = f"{prefix}.xml"

    if not os.path.exists(xml_path):
        return None, "Not Found"

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()

        # Locate calculation type
        calc_type_node = root.find(".//calculation")
        calc_type = calc_type_node.text.strip().lower() if calc_type_node is not None else "unknown"

        # Check exit status
        exit_status_node = root.find(".//exit_status")
        status = "Success" if (exit_status_node is not None and exit_status_node.text == "0") else "Failed"

        return calc_type, status

    except Exception:
        return "error", "Corrupted"
