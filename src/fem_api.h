#pragma once

/**
 * @file fem_api.h
 * @brief C API for FEM solver - enables cross-language interoperability
 *
 * This API provides a C interface to the C++ FEM solver, allowing it to be
 * called from Python, C, or other languages via FFI/ctypes.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Opaque handle to FEM model
typedef void* FEMModelHandle;

/**
 * @brief Create a new FEM model
 * @return Handle to the created model, or NULL on failure
 */
FEMModelHandle FEM_CreateModel();

/**
 * @brief Destroy a FEM model and free all resources
 * @param model Handle to the model
 */
void FEM_DestroyModel(FEMModelHandle model);

/**
 * @brief Add a node to the model
 * @param model Handle to the model
 * @param type Node type (1=TRUSS_NODE, 2=FRAME_NODE)
 * @param x X coordinate
 * @param y Y coordinate
 * @return Node index, or -1 on failure
 */
int FEM_AddNode(FEMModelHandle model, int type, double x, double y);

/**
 * @brief Add an element to the model
 * @param model Handle to the model
 * @param type Element type (1=TRUSS, 2=FRAME)
 * @param node1 First node index
 * @param node2 Second node index
 * @param section Section index
 * @param material Material index
 * @return Element index, or -1 on failure
 */
int FEM_AddElement(FEMModelHandle model, int type, int node1, int node2,
                   int section, int material);

/**
 * @brief Add a material to the model
 * @param model Handle to the model
 * @param E Young's modulus
 * @param mu Poisson's ratio
 * @param alpha Thermal expansion coefficient
 * @return Material index, or -1 on failure
 */
int FEM_AddMaterial(FEMModelHandle model, double E, double mu, double alpha);

/**
 * @brief Add a section to the model
 * @param model Handle to the model
 * @param A Cross-sectional area
 * @param Iz Moment of inertia
 * @param H Section height
 * @return Section index, or -1 on failure
 */
int FEM_AddSection(FEMModelHandle model, double A, double Iz, double H);

/**
 * @brief Add a load to the model
 * @param model Handle to the model
 * @param type Load type
 * @param direction Load direction (0=X, 1=Y, 2=R)
 * @param value Load value
 * @param elem Element index (-1 for node load)
 * @param node Node index (-1 for element load)
 * @param position Position along element
 * @param T0 Temperature at bottom
 * @param T1 Temperature at top
 * @return Load index, or -1 on failure
 */
int FEM_AddLoad(FEMModelHandle model, int type, int direction, double value,
                int elem, int node, double position, double T0, double T1);

/**
 * @brief Add a constraint to the model
 * @param model Handle to the model
 * @param node Node index
 * @param dofX X DOF constraint (-1=constrained, 0=free)
 * @param dofY Y DOF constraint (-1=constrained, 0=free)
 * @param dofR R DOF constraint (-1=constrained, 0=free)
 * @return 0 on success, -1 on failure
 */
int FEM_AddConstraint(FEMModelHandle model, int node, int dofX, int dofY, int dofR);

/**
 * @brief Load model from file
 * @param model Handle to the model
 * @param filename Path to input file
 * @return 0 on success, -1 on failure
 */
int FEM_LoadFromFile(FEMModelHandle model, const char* filename);

/**
 * @brief Save model to file
 * @param model Handle to the model
 * @param filename Path to output file
 * @return 0 on success, -1 on failure
 */
int FEM_SaveToFile(FEMModelHandle model, const char* filename);

/**
 * @brief Solve the FEM problem
 * @param model Handle to the model
 * @return 0 on success, -1 on failure
 */
int FEM_Solve(FEMModelHandle model);

/**
 * @brief Export results to file
 * @param model Handle to the model
 * @param filename Path to results file
 * @return 0 on success, -1 on failure
 */
int FEM_ExportResults(FEMModelHandle model, const char* filename);

/**
 * @brief Get node displacements
 * @param model Handle to the model
 * @param count Output: number of values (3 * number of nodes)
 * @return Pointer to displacement array [Ux0, Uy0, Rz0, Ux1, ...], or NULL
 * @note The returned pointer is valid until the model is destroyed or solved again
 */
const double* FEM_GetDisplacements(FEMModelHandle model, int* count);

/**
 * @brief Get element end forces
 * @param model Handle to the model
 * @param count Output: number of values (6 * number of elements)
 * @return Pointer to forces array [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j, ...], or NULL
 * @note The returned pointer is valid until the model is destroyed or solved again
 */
const double* FEM_GetForces(FEMModelHandle model, int* count);

/**
 * @brief Get support reactions
 * @param model Handle to the model
 * @param count Output: number of values (3 * number of constrained nodes)
 * @return Pointer to reactions array [Rx0, Ry0, Rz0, Rx1, ...], or NULL
 * @note The returned pointer is valid until the model is destroyed or solved again
 */
const double* FEM_GetReactions(FEMModelHandle model, int* count);

/**
 * @brief Get number of nodes in the model
 * @param model Handle to the model
 * @return Number of nodes, or -1 on failure
 */
int FEM_GetNodeCount(FEMModelHandle model);

/**
 * @brief Get number of elements in the model
 * @param model Handle to the model
 * @return Number of elements, or -1 on failure
 */
int FEM_GetElementCount(FEMModelHandle model);

/**
 * @brief Get total number of DOFs
 * @param model Handle to the model
 * @return Total DOF count, or -1 on failure
 */
int FEM_GetDOFCount(FEMModelHandle model);

/**
 * @brief Get number of free DOFs
 * @param model Handle to the model
 * @return Free DOF count, or -1 on failure
 */
int FEM_GetFreeDOFCount(FEMModelHandle model);

/**
 * @brief Get error message from last operation
 * @return Error message string, or empty string if no error
 */
const char* FEM_GetLastError();

#ifdef __cplusplus
}
#endif
