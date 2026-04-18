/**
 * @file pyfem_bindings.cpp
 * @brief Python bindings for FEM library using pybind11
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "src/fem_core.h"
#include "src/fem_types.h"

namespace py = pybind11;

PYBIND11_MODULE(pyfem, m) {
    m.doc() = "Python bindings for FEM_2d finite element analysis library";

    // Constants
    m.attr("TRUSS") = TRUSS;
    m.attr("FRAME") = FRAME;
    m.attr("TRUSS_NODE") = TRUSS_NODE;
    m.attr("FRAME_NODE") = FRAME_NODE;
    m.attr("FORCE_ON_NODE") = FORCE_ON_NODE;
    m.attr("DIRECT_X") = DIRECT_X;
    m.attr("DIRECT_Y") = DIRECT_Y;
    m.attr("DIRECT_R") = DIRECT_R;

    // Material struct
    py::class_<FEM::Material>(m, "Material")
        .def(py::init<>())
        .def_readwrite("E", &FEM::Material::dE, "Young's modulus")
        .def_readwrite("mu", &FEM::Material::dMu, "Poisson's ratio")
        .def_readwrite("alpha", &FEM::Material::dAlpha, "Thermal expansion coefficient");

    // Section struct
    py::class_<FEM::Section>(m, "Section")
        .def(py::init<>())
        .def_readwrite("A", &FEM::Section::dA, "Cross-sectional area")
        .def_readwrite("Iz", &FEM::Section::dIz, "Moment of inertia")
        .def_readwrite("H", &FEM::Section::dH, "Section height");

    // Node struct
    py::class_<FEM::Node>(m, "Node")
        .def(py::init<>())
        .def_readwrite("type", &FEM::Node::iType, "Node type")
        .def_readwrite("x", &FEM::Node::dX, "X coordinate")
        .def_readwrite("y", &FEM::Node::dY, "Y coordinate");

    // Element struct
    py::class_<FEM::Element>(m, "Element")
        .def(py::init<>())
        .def_readwrite("type", &FEM::Element::iType, "Element type")
        .def_property_readonly("nodes",
            [](const FEM::Element& e) {
                return std::vector<int>{e.iaNode[0], e.iaNode[1]};
            }, "Node indices")
        .def_readwrite("section", &FEM::Element::iSection, "Section index")
        .def_readwrite("material", &FEM::Element::iMaterial, "Material index")
        .def_readonly("length", &FEM::Element::dLength, "Element length")
        .def_readonly("sin", &FEM::Element::dSin, "Direction sine")
        .def_readonly("cos", &FEM::Element::dCos, "Direction cosine");

    // Model class
    py::class_<FEM::Model>(m, "Model")
        .def(py::init<>(), "Create a new FEM model")

        // Model building methods
        .def("add_node", &FEM::Model::addNode,
             py::arg("type"), py::arg("x"), py::arg("y"),
             "Add a node to the model\n\n"
             "Args:\n"
             "    type: Node type (TRUSS_NODE=1 or FRAME_NODE=2)\n"
             "    x: X coordinate\n"
             "    y: Y coordinate\n"
             "Returns:\n"
             "    Node index")

        .def("add_element", &FEM::Model::addElement,
             py::arg("type"), py::arg("node1"), py::arg("node2"),
             py::arg("section"), py::arg("material"),
             "Add an element to the model\n\n"
             "Args:\n"
             "    type: Element type (TRUSS=1 or FRAME=2)\n"
             "    node1: First node index\n"
             "    node2: Second node index\n"
             "    section: Section index\n"
             "    material: Material index\n"
             "Returns:\n"
             "    Element index")

        .def("add_material", &FEM::Model::addMaterial,
             py::arg("E"), py::arg("mu"), py::arg("alpha"),
             "Add a material to the model\n\n"
             "Args:\n"
             "    E: Young's modulus (Pa)\n"
             "    mu: Poisson's ratio\n"
             "    alpha: Thermal expansion coefficient\n"
             "Returns:\n"
             "    Material index")

        .def("add_section", &FEM::Model::addSection,
             py::arg("A"), py::arg("Iz"), py::arg("H"),
             "Add a section to the model\n\n"
             "Args:\n"
             "    A: Cross-sectional area (m²)\n"
             "    Iz: Moment of inertia (m⁴)\n"
             "    H: Section height (m)\n"
             "Returns:\n"
             "    Section index")

        .def("add_load", &FEM::Model::addLoad,
             py::arg("type"), py::arg("direction"), py::arg("value"),
             py::arg("elem") = -1, py::arg("node") = -1,
             py::arg("position") = 0.0, py::arg("T0") = 0.0, py::arg("T1") = 0.0,
             "Add a load to the model\n\n"
             "Args:\n"
             "    type: Load type\n"
             "    direction: Load direction (DIRECT_X=0, DIRECT_Y=1, DIRECT_R=2)\n"
             "    value: Load value\n"
             "    elem: Element index (-1 for node load)\n"
             "    node: Node index (-1 for element load)\n"
             "    position: Position along element\n"
             "    T0: Temperature at bottom\n"
             "    T1: Temperature at top\n"
             "Returns:\n"
             "    Load index")

        .def("add_constraint", &FEM::Model::addConstraint,
             py::arg("node"), py::arg("dof_x"), py::arg("dof_y"), py::arg("dof_r"),
             "Add a constraint to the model\n\n"
             "Args:\n"
             "    node: Node index\n"
             "    dof_x: X DOF (-1=constrained, 0=free)\n"
             "    dof_y: Y DOF (-1=constrained, 0=free)\n"
             "    dof_r: R DOF (-1=constrained, 0=free)\n"
             "Returns:\n"
             "    Constraint index")

        // File I/O
        .def("load_from_file", &FEM::Model::loadFromFile,
             py::arg("filename"),
             "Load model from file")

        .def("save_to_file", &FEM::Model::saveToFile,
             py::arg("filename"),
             "Save model to file")

        .def("export_results", &FEM::Model::exportResults,
             py::arg("filename"),
             "Export results to file")

        // Analysis
        .def("solve", &FEM::Model::solve,
             "Solve the FEM problem\n\n"
             "Returns:\n"
             "    True if successful, False otherwise")

        // Results access
        .def("get_displacements", [](FEM::Model& model) {
            int count;
            double* disp = model.getDisplacements(count);
            if (!disp) return py::array_t<double>();
            return py::array_t<double>(count, disp);
        }, "Get node displacements as numpy array")

        .def("get_forces", [](FEM::Model& model) {
            int count;
            double* forces = model.getForces(count);
            if (!forces) return py::array_t<double>();
            return py::array_t<double>(count, forces);
        }, "Get element end forces as numpy array")

        .def("get_reactions", [](FEM::Model& model) {
            int count;
            double* reactions = model.getReactions(count);
            if (!reactions) return py::array_t<double>();
            return py::array_t<double>(count, reactions);
        }, "Get support reactions as numpy array")

        // Getters
        .def_property_readonly("node_count", &FEM::Model::getNodeCount,
                              "Number of nodes in the model")
        .def_property_readonly("element_count", &FEM::Model::getElementCount,
                              "Number of elements in the model")
        .def_property_readonly("dof_count", &FEM::Model::getDOFCount,
                              "Total number of DOFs")
        .def_property_readonly("free_dof_count", &FEM::Model::getFreeDOFCount,
                              "Number of free DOFs");

    // Version info
    m.attr("__version__") = "1.0.0";
}
