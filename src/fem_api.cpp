#include "fem_api.h"
#include "fem_core.h"
#include <string>
#include <cstring>

// Error handling
static std::string g_lastError;

void setLastError(const std::string& error) {
    g_lastError = error;
}

const char* FEM_GetLastError() {
    return g_lastError.c_str();
}

// Model management
FEMModelHandle FEM_CreateModel() {
    try {
        FEM::Model* model = new FEM::Model();
        setLastError("");
        return static_cast<void*>(model);
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to create model: ") + e.what());
        return nullptr;
    }
}

void FEM_DestroyModel(FEMModelHandle model) {
    if (model) {
        delete static_cast<FEM::Model*>(model);
    }
}

// Model building functions
int FEM_AddNode(FEMModelHandle model, int type, double x, double y) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        int idx = m->addNode(type, x, y);
        setLastError("");
        return idx;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to add node: ") + e.what());
        return -1;
    }
}

int FEM_AddElement(FEMModelHandle model, int type, int node1, int node2,
                   int section, int material) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        int idx = m->addElement(type, node1, node2, section, material);
        setLastError("");
        return idx;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to add element: ") + e.what());
        return -1;
    }
}

int FEM_AddMaterial(FEMModelHandle model, double E, double mu, double alpha) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        int idx = m->addMaterial(E, mu, alpha);
        setLastError("");
        return idx;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to add material: ") + e.what());
        return -1;
    }
}

int FEM_AddSection(FEMModelHandle model, double A, double Iz, double H) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        int idx = m->addSection(A, Iz, H);
        setLastError("");
        return idx;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to add section: ") + e.what());
        return -1;
    }
}

int FEM_AddLoad(FEMModelHandle model, int type, int direction, double value,
                int elem, int node, double position, double T0, double T1) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        int idx = m->addLoad(type, direction, value, elem, node, position, T0, T1);
        setLastError("");
        return idx;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to add load: ") + e.what());
        return -1;
    }
}

int FEM_AddConstraint(FEMModelHandle model, int node, int dofX, int dofY, int dofR) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        int idx = m->addConstraint(node, dofX, dofY, dofR);
        setLastError("");
        return idx;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to add constraint: ") + e.what());
        return -1;
    }
}

// File I/O
int FEM_LoadFromFile(FEMModelHandle model, const char* filename) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    if (!filename) {
        setLastError("Invalid filename");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        bool success = m->loadFromFile(filename);
        if (success) {
            setLastError("");
            return 0;
        } else {
            setLastError("Failed to load model from file");
            return -1;
        }
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to load file: ") + e.what());
        return -1;
    }
}

int FEM_SaveToFile(FEMModelHandle model, const char* filename) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    if (!filename) {
        setLastError("Invalid filename");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        bool success = m->saveToFile(filename);
        if (success) {
            setLastError("");
            return 0;
        } else {
            setLastError("Failed to save model to file");
            return -1;
        }
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to save file: ") + e.what());
        return -1;
    }
}

// Analysis
int FEM_Solve(FEMModelHandle model) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        bool success = m->solve();
        if (success) {
            setLastError("");
            return 0;
        } else {
            setLastError("Solver failed");
            return -1;
        }
    } catch (const std::exception& e) {
        setLastError(std::string("Solver error: ") + e.what());
        return -1;
    }
}

int FEM_ExportResults(FEMModelHandle model, const char* filename) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    if (!filename) {
        setLastError("Invalid filename");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        bool success = m->exportResults(filename);
        if (success) {
            setLastError("");
            return 0;
        } else {
            setLastError("Failed to export results");
            return -1;
        }
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to export: ") + e.what());
        return -1;
    }
}

// Results access
const double* FEM_GetDisplacements(FEMModelHandle model, int* count) {
    if (!model || !count) {
        setLastError("Invalid parameters");
        if (count) *count = 0;
        return nullptr;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        double* disp = m->getDisplacements(*count);
        setLastError("");
        return disp;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get displacements: ") + e.what());
        *count = 0;
        return nullptr;
    }
}

const double* FEM_GetForces(FEMModelHandle model, int* count) {
    if (!model || !count) {
        setLastError("Invalid parameters");
        if (count) *count = 0;
        return nullptr;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        double* forces = m->getForces(*count);
        setLastError("");
        return forces;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get forces: ") + e.what());
        *count = 0;
        return nullptr;
    }
}

const double* FEM_GetReactions(FEMModelHandle model, int* count) {
    if (!model || !count) {
        setLastError("Invalid parameters");
        if (count) *count = 0;
        return nullptr;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        double* reactions = m->getReactions(*count);
        setLastError("");
        return reactions;
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get reactions: ") + e.what());
        *count = 0;
        return nullptr;
    }
}

// Model queries
int FEM_GetNodeCount(FEMModelHandle model) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        setLastError("");
        return m->getNodeCount();
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get node count: ") + e.what());
        return -1;
    }
}

int FEM_GetElementCount(FEMModelHandle model) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        setLastError("");
        return m->getElementCount();
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get element count: ") + e.what());
        return -1;
    }
}

int FEM_GetDOFCount(FEMModelHandle model) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        setLastError("");
        return m->getDOFCount();
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get DOF count: ") + e.what());
        return -1;
    }
}

int FEM_GetFreeDOFCount(FEMModelHandle model) {
    if (!model) {
        setLastError("Invalid model handle");
        return -1;
    }
    try {
        FEM::Model* m = static_cast<FEM::Model*>(model);
        setLastError("");
        return m->getFreeDOFCount();
    } catch (const std::exception& e) {
        setLastError(std::string("Failed to get free DOF count: ") + e.what());
        return -1;
    }
}
