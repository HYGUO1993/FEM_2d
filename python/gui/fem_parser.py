import re
import os

class FEMModelData:
    """Contains all data necessary to represent and visualize a FEM model."""
    def __init__(self):
        self.nodes = []           # (type, x, y)
        self.elements = []        # (type, node_i, node_j, section_idx, mat_idx)
        self.constraints = {}     # node_idx -> (cx, cy, cr) where < 0 means constrained
        self.loads = []           # (type, dir, value, elem_idx, node_idx, pos, t0, t1)
        self.materials = []       # (E, mu, alpha)
        self.sections = []        # (A, Iz, H)
        
        # Results data
        self.displacements = {}   # node_idx -> (ux, uy, rz)
        self.end_forces = {}      # elem_idx -> (fxi, fyi, mzi, fxj, fyj, mzj)
        self.reactions = {}       # node_idx -> (rx, ry, rz)
        
        self.is_solved = False

def parse_input_file(filepath: str) -> FEMModelData:
    """Parses a FEM_2d input text file."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")
        
    model = FEMModelData()
    
    with open(filepath, 'r', encoding='utf-8') as f:
        # Filter empty lines and comments
        lines = [line.strip() for line in f.readlines() if line.strip() and not line.strip().startswith('#')]
        
    if not lines:
        raise ValueError("Empty input file")
        
    # Line 0: headers
    headers = list(map(int, lines[0].split()))
    if len(headers) < 6:
        raise ValueError("Invalid header in input file")
        
    n_nodes, n_constraints, n_elements, n_materials, n_sections, n_loads = headers[:6]
    
    idx = 1
    # Nodes
    for _ in range(n_nodes):
        if idx >= len(lines): break
        parts = lines[idx].split()
        if len(parts) >= 3:
            model.nodes.append((int(parts[0]), float(parts[1]), float(parts[2])))
        idx += 1
        
    # Constraints
    for _ in range(n_constraints):
        if idx >= len(lines): break
        parts = lines[idx].split()
        if len(parts) >= 4:
            node_id = int(parts[0])
            model.constraints[node_id] = (int(parts[1]), int(parts[2]), int(parts[3]))
        idx += 1
        
    # Elements
    for _ in range(n_elements):
        if idx >= len(lines): break
        parts = list(map(int, lines[idx].split()))
        if len(parts) >= 4:
            t = parts[0]
            i0 = parts[1]
            i1 = parts[2]
            sec = parts[3]
            mat = parts[4] if len(parts) > 4 else 0
            model.elements.append((t, i0, i1, sec, mat))
        idx += 1
        
    # Materials
    for _ in range(n_materials):
        if idx >= len(lines): break
        parts = list(map(float, lines[idx].split()))
        if len(parts) >= 3:
            model.materials.append((parts[0], parts[1], parts[2]))
        idx += 1
        
    # Sections
    for _ in range(n_sections):
        if idx >= len(lines): break
        parts = list(map(float, lines[idx].split()))
        if len(parts) >= 3:
            model.sections.append((parts[0], parts[1], parts[2]))
        idx += 1
        
    # Loads
    for _ in range(n_loads):
        if idx >= len(lines): break
        parts = lines[idx].split()
        if len(parts) >= 8:
            model.loads.append((
                int(parts[0]), int(parts[1]), float(parts[2]), int(parts[3]), 
                int(parts[4]), float(parts[5]), float(parts[6]), float(parts[7])
            ))
        idx += 1
        
    return model

def parse_results_file(filepath: str, existing_model: FEMModelData = None) -> FEMModelData:
    """Parses the Results.dat file produced by fem_run.exe."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Results file not found: {filepath}")
        
    model = existing_model if existing_model else FEMModelData()
    
    with open(filepath, "r", encoding="utf-8") as f:
        lines = [ln.rstrip("\\n") for ln in f.readlines()]
        
    idx = 0
    n_nodes = 0
    
    # We use raw strings for regex to avoid double escaping problems: r"..."
    while idx < len(lines):
        line = lines[idx].strip()
        
        if line.startswith("Total nodes:"):
            m = re.search(r"Total nodes:\s+(\d+)", line)
            if m:
                n_nodes = int(m.group(1))
                
        # If no existing model, we could parse notes/elements/etc from output (skipped for brevity, 
        # normally we have existing_model from parse_input_file, but occasionally the user just 
        # has the Results.dat or existing_model wasn't provided, so we do basic reading).
        if line.startswith("Node data:") and not existing_model:
            idx += 2
            count = 0
            while idx < len(lines) and count < n_nodes:
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    continue
                parts = t.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    model.nodes.append((int(parts[0]), float(parts[1]), float(parts[2])))
                    count += 1
                idx += 1
            continue
            
        if line.startswith("Element data:") and not existing_model:
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                parts = t.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    model.elements.append((int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])))
                idx += 1
            continue
            
        if line.startswith("Node Displacements:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                # Find all floats, including scientific notation. Crucially, use raw string format and proper negative sign handling
                nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", t)
                if len(nums) >= 4:
                    i = int(float(nums[0]))
                    ux = float(nums[1])
                    uy = float(nums[2])
                    rz = float(nums[3])
                    model.displacements[i] = (ux, uy, rz)
                idx += 1
            continue
            
        if line.startswith("Element End Forces:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", t)
                if len(nums) >= 7:
                    e = int(float(nums[0]))
                    vals = list(map(float, nums[1:7]))
                    model.end_forces[e] = vals
                idx += 1
            continue
            
        if line.startswith("Support Reactions:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", t)
                if len(nums) >= 4:
                    node = int(float(nums[0]))
                    rx = float(nums[1])
                    ry = float(nums[2])
                    rz = float(nums[3])
                    model.reactions[node] = (rx, ry, rz)
                idx += 1
            continue
            
        idx += 1
        
    model.is_solved = bool(model.displacements)
    return model
