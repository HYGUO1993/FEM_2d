#!/usr/bin/env python3
"""
Test script for FEM Python API

This script tests the Python bindings/API by creating a simple cantilever beam
problem and solving it.
"""

import sys
import os

# Add python directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    from fem2d import Model, FRAME_NODE, FRAME, FORCE_ON_NODE, DIRECT_Y
    print(f"✓ Successfully imported fem2d")
    print(f"  Backend: {sys.modules['fem2d']._BACKEND}")
except ImportError as e:
    print(f"✗ Failed to import fem2d: {e}")
    print("  Make sure to build the library first:")
    print("  $ cd build && cmake .. -DBUILD_SHARED_LIB=ON && make")
    sys.exit(1)


def test_cantilever_beam():
    """
    Test case: Cantilever beam with point load at free end

    Geometry:
        Node 0 (0, 0) -------- Node 1 (1, 0)
        Fixed                   Load: -1000N (Y)

    Material: Steel E=210 GPa
    Section: Rectangular 100x100 mm
    """
    print("\n" + "="*60)
    print("Test: Cantilever Beam with Point Load")
    print("="*60)

    # Create model
    model = Model()
    print("✓ Model created")

    # Add material (Steel)
    mat_idx = model.add_material(
        E=210e9,      # Young's modulus (Pa)
        mu=0.3,       # Poisson's ratio
        alpha=1.2e-5  # Thermal expansion coefficient
    )
    print(f"✓ Material added (index: {mat_idx})")

    # Add section (rectangular 100x100 mm)
    A = 0.01  # 100mm * 100mm = 0.01 m²
    Iz = 8.333e-6  # Moment of inertia for rectangular section
    H = 0.1  # Height 100mm
    sec_idx = model.add_section(A=A, Iz=Iz, H=H)
    print(f"✓ Section added (index: {sec_idx})")

    # Add nodes
    n0 = model.add_node(FRAME_NODE, 0.0, 0.0)  # Fixed end
    n1 = model.add_node(FRAME_NODE, 1.0, 0.0)  # Free end
    print(f"✓ Nodes added: {n0}, {n1}")

    # Add element
    e0 = model.add_element(FRAME, n0, n1, sec_idx, mat_idx)
    print(f"✓ Element added (index: {e0})")

    # Add constraint (fixed support at node 0)
    c0 = model.add_constraint(n0, -1, -1, -1)  # All DOFs constrained
    print(f"✓ Constraint added (index: {c0})")

    # Add load (point load at node 1)
    l0 = model.add_load(
        load_type=FORCE_ON_NODE,
        direction=DIRECT_Y,
        value=-1000.0,  # -1000 N downward
        node=n1
    )
    print(f"✓ Load added (index: {l0})")

    # Solve
    print("\nSolving...")
    if model.solve():
        print("✓ Solution successful!")

        # Get results
        try:
            displacements = model.get_displacements()
            if displacements is not None:
                print(f"\nDisplacements (DOF values):")
                print(f"  Total DOFs: {len(displacements)}")
                if len(displacements) >= 6:
                    print(f"  Node 0: Ux={displacements[0]:.6e}, Uy={displacements[1]:.6e}, Rz={displacements[2]:.6e}")
                    print(f"  Node 1: Ux={displacements[3]:.6e}, Uy={displacements[4]:.6e}, Rz={displacements[5]:.6e}")

                    # Theoretical solution for cantilever beam with point load at end:
                    # δ = PL³/(3EI)
                    P = 1000.0
                    L = 1.0
                    E = 210e9
                    I = 8.333e-6
                    delta_theory = (P * L**3) / (3 * E * I)
                    print(f"\n  Theoretical deflection: {delta_theory:.6e} m")
                    print(f"  Computed deflection: {abs(displacements[4]):.6e} m")
                    error = abs(abs(displacements[4]) - delta_theory) / delta_theory * 100
                    print(f"  Relative error: {error:.2f}%")
        except Exception as e:
            print(f"  Note: Could not retrieve displacements: {e}")

        return True
    else:
        print("✗ Solution failed!")
        try:
            error = model.get_last_error()
            if error:
                print(f"  Error: {error}")
        except:
            pass
        return False


def test_file_io():
    """Test loading model from file."""
    print("\n" + "="*60)
    print("Test: Load Model from File")
    print("="*60)

    # Find test file
    test_files = [
        os.path.join(os.path.dirname(__file__), '..', 'test05.txt'),
        os.path.join(os.path.dirname(__file__), '..', 'test_beam.txt'),
    ]

    for test_file in test_files:
        if not os.path.exists(test_file):
            continue

        print(f"\nLoading: {os.path.basename(test_file)}")

        model = Model()
        if model.load_from_file(test_file):
            print(f"✓ Model loaded successfully")
            print(f"  Nodes: {model.node_count}")

            # Solve
            if model.solve():
                print(f"✓ Solution successful")

                # Export results
                result_file = test_file.replace('.txt', '_python_test.dat')
                if model.export_results(result_file):
                    print(f"✓ Results exported to: {os.path.basename(result_file)}")
                    return True
            else:
                print(f"✗ Solution failed")
        else:
            print(f"✗ Failed to load model")

    return False


def main():
    """Run all tests."""
    print("\n" + "╔" + "═"*58 + "╗")
    print("║" + " "*15 + "FEM Python API Tests" + " "*23 + "║")
    print("╚" + "═"*58 + "╝")

    results = []

    # Test 1: Programmatic model creation
    try:
        results.append(("Cantilever Beam Test", test_cantilever_beam()))
    except Exception as e:
        print(f"\n✗ Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Cantilever Beam Test", False))

    # Test 2: File I/O
    try:
        results.append(("File I/O Test", test_file_io()))
    except Exception as e:
        print(f"\n✗ Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("File I/O Test", False))

    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {name}")

    total = len(results)
    passed = sum(1 for _, p in results if p)
    print(f"\nTotal: {passed}/{total} tests passed")

    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
