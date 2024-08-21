import re
import math

def count_functional_groups(smiles):
    # Define the functional groups with their respective SMILES patterns
    functional_groups = {
        'aldehyde': r'C=O(?!C)',  # Aldehyde: C=O not followed by C
        'ketone': r'C\(=O\)C',  # Ketone: C(=O)C
        'carboxylic_acid': r'C\(=O\)\[O-\]',  # Carboxylic Acid: C(=O)[O-]
        'phosphate': r'OP\(=O\)\(\[O-\]\)\[O-\]',  # Phosphate: OP(=O)([O-])[O-]
        'phosphate(c)': r'OP\(=O\)\(\[O-\]\)(?!\[O-\])',  # Consecutive phosphate: OP(=O)([O-]), not followed by [O-]
        'acyl_phosphate': r'C\(=O\)OP\(=O\)\(\[O-\]\)\[O-\]',  # Acyl Phosphate: C(=O)OP(=O)([O-])[O-]
        'alkene': r'C=C', # Alkene: C=C
        'carbons': r'C', # Carbon: C
        'ether': r'(?<!C\(=O\))O1|O3',#|O2???  # Ether: O bonded to two carbons (can be O1 or O3 in SMILES)
        'amine': r'N(?![\w\[\(])',  # Amine: N not followed by any alphanumeric characters, brackets, or parentheses
        'purine': r'1=NC\(=C2C\(=N1\)N\(C=N2\)', # Purine ring
        'lactone': r'\(C\(=O\)O1\)', # Lactone: (C(=O)O1)
        
    }

    # Initialize a dictionary to store counts
    counts = {group: 0 for group in functional_groups}

    # Count occurrences of each functional group using regex
    for group, pattern in functional_groups.items():
        matches = re.finditer(pattern, smiles)
        counts[group] = sum(1 for _ in matches)

    # Adjust phosphate count to exclude acyl phosphates
    counts['phosphate'] -= counts['acyl_phosphate']

    # Count total oxygens
    total_oxygens = len(re.findall(r'O', smiles))

    # Calculate the total number of oxygens contributed by each functional group
    oxygens_in_other_groups = (counts['aldehyde'] +  # Each aldehyde contributes 1 oxygen
                                counts['ketone'] +  # Each ketone contributes 1 oxygen
                                counts['carboxylic_acid'] * 2 +  # Each carboxylic acid contributes 2 oxygens
                                counts['acyl_phosphate'] * 5 +  # Each acyl phosphate contributes 5 oxygens
                                counts['phosphate'] * 4 +  # Each phosphate contributes 4 oxygens
                                counts['phosphate(c)'] * 3 +  # Each consecutive phosphate contributes 3 oxygens
                                counts['ether'] +  # Each ether contributes 1 oxygen
                                counts['lactone'] * 2) # Each lactone contributes 2 oxygens

    # Calculate hydroxyls
    total_hydroxyls = total_oxygens - oxygens_in_other_groups

    # Update counts with the corrected hydroxyl count
    counts['hydroxyl'] = total_hydroxyls

    return counts

def calculate_entropy(counts):
    # Calculate the Shannon entropy based on the counts of functional groups and carbons
    total_groups = sum(counts.values())
    entropy = 0.0

    for count in counts.values():
        if count > 0:
            probability = count / total_groups
            entropy -= probability * math.log2(probability)

    return entropy

# Example usage:
smiles = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N"
counts = count_functional_groups(smiles)
entropy = calculate_entropy(counts)

# Print the results
print("Counts of functional groups and carbons:", counts)
print("Shannon Entropy:", entropy)

