import numpy as np

def print_header(title, width=40):
    print("┌" + "─" * (width-2) + "┐")
    print(f"│ {title.center(width-4)} │")
    print("├" + "─" * (width-2) + "┤")

def print_footer(width=40):
    print("└" + "─" * (width-2) + "┘")

def print_results(title, values, value_label="Value", unit=""):
    """
    Print numerical results in a formatted box with aligned columns.
    
    Args:
        title (str): Header title
        values (list): List of numerical values
        value_label (str): Description of values (e.g., "Complex")
        unit (str): Measurement unit (e.g., "kJ/mol")
    """
    max_width = max(len(f"{v:.2f}") for v in values) if values else 10
    width = max(40, len(title) + 8, len(value_label) + max_width + 12)
    
    print_header(title, width)
    
    if not values:
        print(f"│ {'No results'.center(width-4)} │")
    else:
        for i, val in enumerate(values, 1):
            if isinstance(val, (np.floating, float)):
                formatted_val = f"{val:.2f}" if abs(val) >= 0.01 else f"{val:.2e}"
            else:
                formatted_val = str(val)
                
            line = f"│ {value_label} {i}: {formatted_val:>{max_width}} {unit} │"
            print(line.ljust(width-1) + "│")
    
    print_footer(width)

def print_equilibrium(result):
    """Specialized formatter for equilibrium results"""
    width = 40
    # Create lines with aligned columns
    lines = [
        f"Free Ligand: {result['free_ligand']:.3f} M",
        f"Free Metal: {result['free_metal']:.2e} M"
    ]
    
    # Add complexes with consistent numbering
    for i, conc in enumerate(result['complexes'], 1):
        if conc < 0.01:
            lines.append(f"Complex {i}: {conc:.2e} M")
        else:
            lines.append(f"Complex {i}: {conc:.3f} M")
    
    # Find longest label length
    max_label_len = max(len(line.split(":")[0]) for line in lines)
    
    # Build formatted output
    print("┌─────────────────────────────────┐")
    print("│      Equilibrium Concentrations │")
    print("├─────────────────────────────────┤")
    
    for line in lines:
        label, value = line.split(":")
        padded_label = label.ljust(max_label_len + 2)
        padded_line = f"│ {padded_label}: {value.strip():<16} │"
        print(padded_line)
    
    print("└─────────────────────────────────┘")