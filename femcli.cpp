// femcli.cpp — FEM_2d 命令行工具(JSON 模式 + txt 兼容模式)
//
// 用法:
//   femcli solve model.json -o result.json          # JSON 输入 → JSON 输出(推荐)
//   femcli solve -i model.txt -o result.json --legacy-txt   # txt 输入 → JSON 输出
//   femcli validate model.json                      # 仅校验模型,不求解
//
// 设计目标:统一模型格式 JSON,让 GUI / CLI / LLM 三个入口共享同一份数据。

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "barsystem.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;
using namespace std;

// ---------------------------------------------------------------------------
// 枚举 ↔ 字符串映射
// ---------------------------------------------------------------------------
static int NodeTypeFromStr(const string& s)
{
	if (s == "truss") return TRUSS_NODE;
	return FRAME_NODE;
}

static int ElemTypeFromStr(const string& s)
{
	if (s == "truss") return TRUSS;
	return FRAME;
}

static int LoadTypeFromStr(const string& s)
{
	if (s == "nodalForce")               return FORCE_ON_NODE;
	if (s == "lateralForce")             return LATERAL_FORCE;
	if (s == "lateralUniformPressure")   return LATERAL_UNIFORM_PRESSURE;
	if (s == "momentOnPoint")            return MOMENT_ON_A_POINT;
	if (s == "lateralLinearlyPressure")  return LATERAL_LINEARLY_PRESSURE;
	if (s == "axialPressure")            return AXIAL_PRESSURE;
	if (s == "axialForce")               return AXIAL_FORCE;
	if (s == "momentOnBeam")             return MOMENT_ON_BEAM;
	if (s == "temperature")              return TEMPERATURE;
	if (s == "supportMove")              return SUPPORT_MOVE;
	return -1; // 未知
}

static int DirectFromStr(const string& s)
{
	if (s == "x")  return DIRECT_X;
	if (s == "y")  return DIRECT_Y;
	if (s == "rz") return DIRECT_R;
	return 0;
}

// ---------------------------------------------------------------------------
// JSON → 内部结构
// ---------------------------------------------------------------------------
static bool LoadJsonModel(const json& j, std::string& err,
                          int& nNode, int& nCons, int& nElem,
                          int& nMat, int& nSect, int& nLoad,
                          vector<Node>& nodes, vector<ConstrainedNode>& cons,
                          vector<Element>& elems, vector<Material>& mats,
                          vector<Section>& sects, vector<Load>& loads)
{
	// --- 节点 ---
	if (!j.contains("nodes") || !j["nodes"].is_array() || j["nodes"].empty()) {
		err = "缺少 nodes 数组"; return false;
	}
	const auto& jn = j["nodes"];
	nNode = (int)jn.size();
	nodes.resize(nNode);
	for (int i = 0; i < nNode; ++i) {
		const auto& n = jn[i];
		nodes[i].iType = NodeTypeFromStr(n.value("type", "frame"));
		nodes[i].dX = n.value("x", 0.0);
		nodes[i].dY = n.value("y", 0.0);
	}

	// --- 约束 ---
	nCons = 0;
	cons.clear();
	if (j.contains("constraints") && j["constraints"].is_array()) {
		nCons = (int)j["constraints"].size();
		cons.resize(nCons);
		for (int i = 0; i < nCons; ++i) {
			const auto& c = j["constraints"][i];
			cons[i].iNode = c.value("node", 0);
			cons[i].iaConstrainedDOF[0] = 0;
			cons[i].iaConstrainedDOF[1] = 0;
			cons[i].iaConstrainedDOF[2] = 0;
			if (c.contains("dofs") && c["dofs"].is_array()) {
				for (const auto& d : c["dofs"]) {
					const string s = d.get<string>();
					if (s == "ux") cons[i].iaConstrainedDOF[0] = -1;
					else if (s == "uy") cons[i].iaConstrainedDOF[1] = -1;
					else if (s == "rz") cons[i].iaConstrainedDOF[2] = -1;
				}
			} else {
				// 兼容旧字段: dx/dy/dr 数字
				cons[i].iaConstrainedDOF[0] = c.value("dx", 0);
				cons[i].iaConstrainedDOF[1] = c.value("dy", 0);
				cons[i].iaConstrainedDOF[2] = c.value("dr", 0);
			}
		}
	}

	// --- 单元 ---
	if (!j.contains("elements") || !j["elements"].is_array() || j["elements"].empty()) {
		err = "缺少 elements 数组"; return false;
	}
	const auto& je = j["elements"];
	nElem = (int)je.size();
	elems.resize(nElem);
	for (int i = 0; i < nElem; ++i) {
		const auto& e = je[i];
		elems[i].iType = ElemTypeFromStr(e.value("type", "frame"));
		elems[i].iaNode[0] = e.value("nodeI", 0);
		elems[i].iaNode[1] = e.value("nodeJ", 0);
		elems[i].iSection = e.value("section", 0);
		elems[i].iMaterial = e.value("material", 0);
	}

	// --- 材料 ---
	if (!j.contains("materials") || !j["materials"].is_array() || j["materials"].empty()) {
		err = "缺少 materials 数组"; return false;
	}
	const auto& jm = j["materials"];
	nMat = (int)jm.size();
	mats.resize(nMat);
	for (int i = 0; i < nMat; ++i) {
		const auto& m = jm[i];
		mats[i].dE = m.value("E", 0.0);
		mats[i].dMu = m.value("mu", 0.0);
		mats[i].dAlpha = m.value("alpha", 0.0);
	}

	// --- 截面 ---
	if (!j.contains("sections") || !j["sections"].is_array() || j["sections"].empty()) {
		err = "缺少 sections 数组"; return false;
	}
	const auto& js = j["sections"];
	nSect = (int)js.size();
	sects.resize(nSect);
	for (int i = 0; i < nSect; ++i) {
		const auto& s = js[i];
		sects[i].dA = s.value("A", 0.0);
		sects[i].dIz = s.value("Iz", 0.0);
		sects[i].dH = s.value("height", 0.0);
	}

	// --- 荷载 ---
	nLoad = 0;
	loads.clear();
	if (j.contains("loads") && j["loads"].is_array()) {
		nLoad = (int)j["loads"].size();
		loads.resize(nLoad);
		for (int i = 0; i < nLoad; ++i) {
			const auto& L = j["loads"][i];
			loads[i].iType = LoadTypeFromStr(L.value("type", "nodalForce"));
			loads[i].iDirect = DirectFromStr(L.value("direction", "y"));
			loads[i].dValue = L.value("value", 0.0);
			loads[i].iLoadedElem = L.value("element", -1);
			loads[i].iLoadedNode = L.value("node", -1);
			loads[i].dPosition = L.value("position", 0.0);
			loads[i].dT0 = L.value("T0", 0.0);
			loads[i].dT1 = L.value("T1", 0.0);
		}
	}
	return true;
}

// ---------------------------------------------------------------------------
// txt → 内部结构(与 fem_run.exe 完全一致,供 --legacy-txt 使用)
// ---------------------------------------------------------------------------
static bool LoadTxtModel(const string& path, std::string& err,
                         int& nNode, int& nCons, int& nElem,
                         int& nMat, int& nSect, int& nLoad,
                         vector<Node>& nodes, vector<ConstrainedNode>& cons,
                         vector<Element>& elems, vector<Material>& mats,
                         vector<Section>& sects, vector<Load>& loads)
{
	ifstream fin(path.c_str());
	if (!fin) { err = "无法打开输入文件: " + path; return false; }
	fin >> nNode >> nCons >> nElem >> nMat >> nSect >> nLoad;
	if (fin.fail() || nNode <= 0 || nElem <= 0 || nMat <= 0 || nSect <= 0 || nLoad < 0) {
		err = "无效的总控数据"; return false;
	}
	nodes.resize(nNode);
	for (int i = 0; i < nNode; ++i)
		fin >> nodes[i].iType >> nodes[i].dX >> nodes[i].dY;
	cons.resize(nCons);
	for (int i = 0; i < nCons; ++i)
		fin >> cons[i].iNode >> cons[i].iaConstrainedDOF[0]
		    >> cons[i].iaConstrainedDOF[1] >> cons[i].iaConstrainedDOF[2];
	elems.resize(nElem);
	for (int i = 0; i < nElem; ++i) {
		fin >> elems[i].iType >> elems[i].iaNode[0] >> elems[i].iaNode[1] >> elems[i].iSection;
		string rest;
		getline(fin, rest);
		istringstream iss(rest);
		int matIndex = 0;
		if (!(iss >> matIndex)) matIndex = 0;
		elems[i].iMaterial = matIndex;
	}
	mats.resize(nMat);
	for (int i = 0; i < nMat; ++i)
		fin >> mats[i].dE >> mats[i].dMu >> mats[i].dAlpha;
	sects.resize(nSect);
	for (int i = 0; i < nSect; ++i)
		fin >> sects[i].dA >> sects[i].dIz >> sects[i].dH;
	loads.resize(nLoad);
	for (int i = 0; i < nLoad; ++i)
		fin >> loads[i].iType >> loads[i].iDirect >> loads[i].dValue
		    >> loads[i].iLoadedElem >> loads[i].iLoadedNode >> loads[i].dPosition
		    >> loads[i].dT0 >> loads[i].dT1;
	return true;
}

// ---------------------------------------------------------------------------
// 内部结构 → JSON 结果
// ---------------------------------------------------------------------------
static json BuildResultJson(int nTotalNode, int nCons, int nElem,
                            int nMat, int nSect, int nLoad,
                            const vector<Node>& nodes,
                            const vector<ConstrainedNode>& cons,
                            const vector<Element>& elems,
                            const vector<Material>& mats,
                            const vector<Section>& sects,
                            const vector<Load>& loads,
                            const double* pDisp, const double* pLoadVect,
                            int nTotalDOF, int nFreeDOF)
{
	json out;
	out["schemaVersion"] = "1.0";
	out["solver"] = "builtin";
	out["status"] = "ok";
	out["message"] = "求解完成";
	out["stats"] = {
		{"nodeCount", nTotalNode},
		{"elementCount", nElem},
		{"totalDOF", nTotalDOF},
		{"freeDOF", nFreeDOF}
	};

	// 节点位移
	json disps = json::array();
	for (int i = 0; i < nTotalNode; ++i) {
		json d;
		d["node"] = i;
		d["ux"] = pDisp[nodes[i].iaDOFIndex[0]];
		d["uy"] = pDisp[nodes[i].iaDOFIndex[1]];
		d["rz"] = pDisp[nodes[i].iaDOFIndex[2]];
		disps.push_back(d);
	}
	out["displacements"] = disps;

	// 单元端力(局部坐标 N V M)
	json efs = json::array();
	for (int i = 0; i < nElem; ++i) {
		json e;
		e["element"] = i;
		e["type"] = (elems[i].iType == TRUSS) ? "truss" : "frame";
		e["nodeI"] = { {"N", elems[i].daEndInterForce[0]},
		               {"V", elems[i].daEndInterForce[1]},
		               {"M", elems[i].daEndInterForce[2]} };
		e["nodeJ"] = { {"N", elems[i].daEndInterForce[3]},
		               {"V", elems[i].daEndInterForce[4]},
		               {"M", elems[i].daEndInterForce[5]} };
		efs.push_back(e);
	}
	out["endForces"] = efs;

	// 支座反力: 约束 DOF 上的 pLoadVect(位移分量叠加后为总荷载=反力)
	json reacs = json::array();
	for (int i = 0; i < nCons; ++i) {
		json r;
		r["node"] = cons[i].iNode;
		r["ux"] = (cons[i].iaConstrainedDOF[0] < 0) ? pLoadVect[nodes[cons[i].iNode].iaDOFIndex[0]] : 0.0;
		r["uy"] = (cons[i].iaConstrainedDOF[1] < 0) ? pLoadVect[nodes[cons[i].iNode].iaDOFIndex[1]] : 0.0;
		r["rz"] = (cons[i].iaConstrainedDOF[2] < 0) ? pLoadVect[nodes[cons[i].iNode].iaDOFIndex[2]] : 0.0;
		reacs.push_back(r);
	}
	out["reactions"] = reacs;
	return out;
}

// ---------------------------------------------------------------------------
// 主流程
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	string cmd, inputPath, outputPath = "result.json";
	bool legacyTxt = false;

	for (int i = 1; i < argc; ++i) {
		string a = argv[i];
		if (a == "solve") cmd = "solve";
		else if (a == "validate") cmd = "validate";
		else if (a == "-i" || a == "--input")  { if (i+1<argc) inputPath = argv[++i]; }
		else if (a == "-o" || a == "--output") { if (i+1<argc) outputPath = argv[++i]; }
		else if (a == "--legacy-txt") legacyTxt = true;
		else if (inputPath.empty() && a[0] != '-') inputPath = a;
	}

	if (cmd.empty() || inputPath.empty()) {
		cerr << "用法: femcli solve <model.json> -o <result.json> [--legacy-txt]" << endl;
		cerr << "      femcli validate <model.json>" << endl;
		return 1;
	}

	int nNode=0, nCons=0, nElem=0, nMat=0, nSect=0, nLoad=0;
	vector<Node> nodes;
	vector<ConstrainedNode> cons;
	vector<Element> elems;
	vector<Material> mats;
	vector<Section> sects;
	vector<Load> loads;
	string err;

	if (legacyTxt) {
		if (!LoadTxtModel(inputPath, err, nNode, nCons, nElem, nMat, nSect, nLoad,
		                  nodes, cons, elems, mats, sects, loads)) {
			cerr << "ERROR: " << err << endl;
			return 1;
		}
	} else {
		ifstream fin(inputPath.c_str());
		if (!fin) { cerr << "ERROR: 无法打开 " << inputPath << endl; return 1; }
		json j;
		try {
			fin >> j;
		} catch (const json::parse_error& pe) {
			cerr << "ERROR: JSON 解析失败: " << pe.what() << endl;
			return 1;
		}
		if (!LoadJsonModel(j, err, nNode, nCons, nElem, nMat, nSect, nLoad,
		                   nodes, cons, elems, mats, sects, loads)) {
			cerr << "ERROR: " << err << endl;
			return 1;
		}
	}

	if (cmd == "validate") {
		cout << "模型校验通过: " << nNode << " 节点, " << nElem << " 单元, "
		     << nLoad << " 荷载" << endl;
		return 0;
	}

	// 求解
	int nFreeDOF = 0;
	int nTotalDOF = DOFIndexCalcu(nFreeDOF, nNode, nCons,
	                              cons.data(), nodes.data());
	vector<double> pDisp(nTotalDOF, 0.0);
	vector<double> pLoadVect(nTotalDOF, 0.0);
	string solveErr;
	bool ok = FemSolveModel(nNode, nCons, nElem, nMat, nSect, nLoad,
	                        nodes.data(), cons.data(), elems.data(),
	                        mats.data(), sects.data(), loads.data(),
	                        pDisp.data(), pLoadVect.data(),
	                        "", false, true, solveErr);
	if (!ok) {
		json errOut;
		errOut["schemaVersion"] = "1.0";
		errOut["solver"] = "builtin";
		errOut["status"] = "error";
		errOut["message"] = solveErr;
		ofstream fout(outputPath.c_str());
		fout << errOut.dump(2) << endl;
		fout.close();
		cerr << "ERROR: " << solveErr << endl;
		return 2;
	}

	json result = BuildResultJson(nNode, nCons, nElem, nMat, nSect, nLoad,
	                              nodes, cons, elems, mats, sects, loads,
	                              pDisp.data(), pLoadVect.data(),
	                              nTotalDOF, nFreeDOF);
	ofstream fout(outputPath.c_str());
	if (!fout) { cerr << "ERROR: 无法写入 " << outputPath << endl; return 1; }
	fout << result.dump(2) << endl;
	fout.close();

	if (argc >= 1 && getenv("FEMCLI_VERBOSE")) {
		cout << "求解完成: " << nNode << " 节点, " << nElem << " 单元, "
		     << nFreeDOF << "/" << nTotalDOF << " DOF" << endl;
	}
	return 0;
}
