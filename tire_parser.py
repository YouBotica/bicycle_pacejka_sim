"""
Parses Pacejka coefficients from a .tir file string content.
"""
import re

def parse_tir_coeffs(tir_content):
    """
    Parses Pacejka coefficients from the string content of a .tir file.

    Args:
        tir_content (str): The full string content of the .tir file.

    Returns:
        dict: A dictionary containing the extracted coefficients, structured by
              section (LONGITUDINAL, LATERAL, ALIGNING). Returns None if
              parsing fails or essential sections are missing.
    """
    coeffs = {
        'LONGITUDINAL': {},
        'LATERAL': {},
        'ALIGNING': {}
        # Add other sections like 'OVERTURNING', 'ROLLING' if needed later
    }
    current_section = None
    section_mapping = {
        '[LONGITUDINAL_COEFFICIENTS]': 'LONGITUDINAL',
        '[LATERAL_COEFFICIENTS]': 'LATERAL',
        '[ALIGNING_COEFFICIENTS]': 'ALIGNING',
        # Add mappings for other sections if needed
    }

    lines = tir_content.splitlines()

    for line in lines:
        line = line.strip()
        if not line or line.startswith('$') or line.startswith('!'):
            continue

        # Check for section headers
        is_section_header = False
        for header, section_key in section_mapping.items():
            if line.upper() == header:
                current_section = section_key
                is_section_header = True
                break
        if is_section_header:
            continue

        # If inside a relevant section, parse key-value pairs
        if current_section:
            # Use regex to handle potential comments after the value
            match = re.match(r'^\s*([A-Za-z0-9_]+)\s*=\s*([-\d\.eE+]+)', line)
            if match:
                key = match.group(1).strip()
                try:
                    value = float(match.group(2).strip())
                    coeffs[current_section][key] = value
                except ValueError:
                    print(f"Warning: Could not parse value for {key} in section {current_section}: {match.group(2)}")
            # Stop parsing when the next section starts
            elif line.startswith('['):
                 current_section = None # Reset section if a new, unmapped one starts

    # Basic validation: Check if essential coefficients were found
    if not coeffs['LONGITUDINAL'] or not coeffs['LATERAL'] or not coeffs['ALIGNING']:
        print("Error: Failed to parse essential coefficient sections.")
        return None

    return coeffs

# Example usage (assuming tir_file_content holds the string from read_file)
# tir_file_content = """... (content of lf_tire.tir) ..."""
# parsed_coefficients = parse_tir_coeffs(tir_file_content)
# if parsed_coefficients:
#     print("Successfully parsed coefficients.")
#     # print(parsed_coefficients['LATERAL']['PCY1']) # Example access
# else:
#     print("Failed to parse coefficients.")