#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "PNI_Calculate.hpp"

using namespace std;
using namespace Eigen;
namespace fs = std::filesystem;

struct Atom {
    string atomType;
    string atomName;
    string resName;
    char chainID;
    int resID;
    double x, y, z;
};

using Residue = vector<Atom>;
using Chain = map<int, Residue>;      // Residue ID -> Residue
using Structure = map<char, Chain>;   // Chain ID -> Chain

string trim(const string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    size_t last = str.find_last_not_of(" \t\r\n");
    if (first == string::npos || last == string::npos) {
        return "";
    }
    return str.substr(first, last - first + 1);
}

Structure ParsePDB(const string& filename) {
    Structure structure;
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Failed to open PDB file: " << filename << endl;
        return structure;
    }

    string line;
    while (getline(file, line)) {
        if (line.size() < 54) {
            continue;
        }
        if (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0) {
            try {
                Atom atom;
                atom.atomType = trim(line.substr(0, 6));
                atom.atomName = trim(line.substr(12, 4));
                atom.resName = trim(line.substr(17, 3));
                atom.chainID = line[21];
                atom.resID = stoi(line.substr(22, 4));
                atom.x = stod(line.substr(30, 8));
                atom.y = stod(line.substr(38, 8));
                atom.z = stod(line.substr(46, 8));
                structure[atom.chainID][atom.resID].push_back(atom);
            } catch (const exception&) {
                continue;
            }
        }
    }
    return structure;
}

map<string, Residue> baseResDB;

fs::path ResolveTemplatePath(const char* argv0) {
    vector<fs::path> candidates = {
        fs::path("data") / "nucleicacid.dat",
        fs::path("nucleicacid.dat")
    };

    if (argv0) {
        error_code ec;
        fs::path exeDir = fs::weakly_canonical(fs::path(argv0), ec).parent_path();
        if (!exeDir.empty()) {
            candidates.push_back(exeDir / "data" / "nucleicacid.dat");
            candidates.push_back(exeDir / "nucleicacid.dat");
            candidates.push_back(exeDir.parent_path() / "data" / "nucleicacid.dat");
        }
    }

    for (const auto& candidate : candidates) {
        if (fs::exists(candidate)) {
            return candidate;
        }
    }
    return fs::path("data") / "nucleicacid.dat";
}

bool LoadBaseResidue(const fs::path& dataPath) {
    ifstream file(dataPath);
    if (!file) {
        cerr << "Failed to open base template file: " << dataPath << endl;
        return false;
    }

    string line;
    int lineCount = 0;

    while (getline(file, line)) {
        lineCount++;
        if (trim(line).empty()) {
            continue;
        }
        istringstream line_ss(line);
        Atom atomParam{};
        if (line_ss >> atomParam.atomName >> atomParam.resName >> atomParam.x >> atomParam.y >> atomParam.z) {
            baseResDB[atomParam.resName].push_back(atomParam);
        } else {
            cerr << "Invalid line in " << dataPath.filename() << ": " << lineCount << "." << endl;
        }
    }

    if (baseResDB.empty()) {
        cerr << "Error: No base residues loaded from " << dataPath << endl;
        return false;
    }

    return true;
}

void HelixFromTmpl(Structure& structure, const string& chainID, int resID_begin, int resID_end, const string& seq) {
    string fullSeq;
    fullSeq.reserve(seq.size());
    for (char c : seq) {
        fullSeq += static_cast<char>(toupper(static_cast<unsigned char>(c)));
    }

    if (chainID.empty() || structure.find(chainID[0]) == structure.end()) {
        cerr << "Error: Chain '" << chainID << "' not found in structure." << endl;
        return;
    }

    const vector<string> backboneAtoms = {
        "P", "OP1", "OP2", "O2'", "O3'", "O5'", "C1'", "C2'", "C3'", "C4'", "C5'", "O4'"
    };

    auto& chain = structure[chainID[0]];

    bool isDNA = false;
    for (const auto& [rID, res] : chain) {
        if (!res.empty() && res[0].resName.find('D') != string::npos) {
            isDNA = true;
            break;
        }
    }

    for (size_t i = 0; i < fullSeq.size(); ++i) {
        int resID = resID_begin + static_cast<int>(i);
        if (resID > resID_end) {
            break;
        }

        auto resIt = chain.find(resID);
        if (resIt == chain.end() || resIt->second.empty()) {
            cerr << "Warning: Residue " << chainID << resID << " not found in structure, skipping." << endl;
            continue;
        }
        auto& residue = resIt->second;

        string targetBase(1, fullSeq[i]);
        auto tmplIt = baseResDB.find(targetBase);
        if (tmplIt == baseResDB.end() || tmplIt->second.empty()) {
            cerr << "Error: Base template not found for '" << targetBase << "', skipping residue " << resID << "." << endl;
            continue;
        }
        const Residue& baseTmpl = tmplIt->second;

        // Preserve sugar-phosphate backbone atoms
        Residue newResidue;
        string newResName = isDNA ? ("D" + targetBase) : targetBase;

        for (const auto& atom : residue) {
            if (find(backboneAtoms.begin(), backboneAtoms.end(), atom.atomName) != backboneAtoms.end()) {
                Atom copiedAtom = atom;
                copiedAtom.resName = newResName;
                copiedAtom.atomType = "ATOM";
                newResidue.push_back(copiedAtom);
            }
        }

        // Compute original base plane normal vector n1
        vector<Atom> origBaseAtoms;
        bool isOrigPurine = (residue[0].resName == "DA" || residue[0].resName == "DG" ||
                             residue[0].resName == "A" || residue[0].resName == "G");
        vector<string> origBaseAtomNames = isOrigPurine
                                               ? vector<string>{"N9", "C4", "C5"}
                                               : vector<string>{"N1", "C2", "N3"};
        for (const string& atomName : origBaseAtomNames) {
            for (const auto& atom : residue) {
                if (atom.atomName == atomName) {
                    origBaseAtoms.push_back(atom);
                    break;
                }
            }
        }
        if (origBaseAtoms.size() < 3) {
            cerr << "Warning: Not enough base atoms to define plane for residue " << resID << ", skipping." << endl;
            continue;
        }

        Vector3d v1(origBaseAtoms[1].x - origBaseAtoms[0].x,
                    origBaseAtoms[1].y - origBaseAtoms[0].y,
                    origBaseAtoms[1].z - origBaseAtoms[0].z);
        Vector3d v2(origBaseAtoms[2].x - origBaseAtoms[0].x,
                    origBaseAtoms[2].y - origBaseAtoms[0].y,
                    origBaseAtoms[2].z - origBaseAtoms[0].z);
        Vector3d n1 = v1.cross(v2).normalized();

        // Compute template base plane normal vector n2 and rotation R1 aligning n2 -> n1
        vector<Atom> newBaseAtoms;
        bool isTmplPurine = (baseTmpl[0].resName == "A" || baseTmpl[0].resName == "G");
        vector<string> newBaseAtomNames = isTmplPurine
                                              ? vector<string>{"N9", "C4", "C5"}
                                              : vector<string>{"N1", "C2", "N3"};
        for (const string& atomName : newBaseAtomNames) {
            for (const auto& atom : baseTmpl) {
                if (atom.atomName == atomName) {
                    newBaseAtoms.push_back(atom);
                    break;
                }
            }
        }
        if (newBaseAtoms.size() < 3) {
            cerr << "Warning: Not enough base atoms in template for " << targetBase << ", skipping." << endl;
            continue;
        }

        Vector3d nv1(newBaseAtoms[1].x - newBaseAtoms[0].x,
                     newBaseAtoms[1].y - newBaseAtoms[0].y,
                     newBaseAtoms[1].z - newBaseAtoms[0].z);
        Vector3d nv2(newBaseAtoms[2].x - newBaseAtoms[0].x,
                     newBaseAtoms[2].y - newBaseAtoms[0].y,
                     newBaseAtoms[2].z - newBaseAtoms[0].z);
        Vector3d n2 = nv1.cross(nv2).normalized();
        Matrix3d R1 = Quaterniond::FromTwoVectors(n2, n1).toRotationMatrix();

        // Check for non-standard HETATM residue; prompt for glycosidic bond atoms if needed
        bool isHet = false;
        for (const Atom& atom : residue) {
            if (atom.atomType == "HETATM") {
                isHet = true;
                break;
            }
        }

        string atom1 = "C1'";
        string atom2 = isOrigPurine ? "N9" : "N1";
        if (isHet) {
            cout << "Residue " << chainID << resID
                 << " is HETATM. Enter two atom names for glycosidic bond alignment:" << endl;
            cout << "First atom (on sugar): ";
            cin >> atom1;
            cout << "Second atom (on base): ";
            cin >> atom2;
        }

        Atom* origC1 = nullptr;
        Atom* origGlycosidicAtom = nullptr;
        for (auto& atom : residue) {
            if (atom.atomName == atom1) {
                origC1 = &atom;
                break;
            }
        }
        for (auto& atom : residue) {
            if (atom.atomName == atom2) {
                origGlycosidicAtom = &atom;
                break;
            }
        }
        if (!origC1 || !origGlycosidicAtom) {
            cerr << "Warning: Missing alignment atoms (" << atom1 << ", " << atom2
                 << ") for residue " << chainID << resID << ", skipping." << endl;
            continue;
        }

        Vector3d oldC1(origC1->x, origC1->y, origC1->z);
        Vector3d oldGlycosidicAtom(origGlycosidicAtom->x, origGlycosidicAtom->y, origGlycosidicAtom->z);

        // Locate template glycosidic bond atoms
        string glycosidicAtomTmpl = isTmplPurine ? "N9" : "N1";
        bool foundTmplC1 = false;
        bool foundTmplGly = false;
        Vector3d newC1 = Vector3d::Zero();
        Vector3d newGlycosidicAtom = Vector3d::Zero();
        for (const auto& atom : baseTmpl) {
            if (atom.atomName == "C1'") {
                newC1 = Vector3d(atom.x, atom.y, atom.z);
                foundTmplC1 = true;
            } else if (atom.atomName == glycosidicAtomTmpl) {
                newGlycosidicAtom = Vector3d(atom.x, atom.y, atom.z);
                foundTmplGly = true;
            }
        }
        if (!foundTmplC1 || !foundTmplGly) {
            cerr << "Warning: Missing glycosidic atoms in template for " << targetBase << ", skipping." << endl;
            continue;
        }

        // Compute in-plane rotation R2 to align glycosidic bond vectors
        Vector3d oldGlycosidicVector = (oldGlycosidicAtom - oldC1).normalized();
        Vector3d newGlycosidicVector = (newGlycosidicAtom - newC1).normalized();
        Vector3d newGlyVecRotated = R1 * newGlycosidicVector;

        auto ProjectOntoPlane = [](const Vector3d& v, const Vector3d& normal) {
            return (v - v.dot(normal) * normal).normalized();
        };
        Vector3d newGlycosidicProj = ProjectOntoPlane(newGlyVecRotated, n1);
        Vector3d oldGlycosidicProj = ProjectOntoPlane(oldGlycosidicVector, n1);

        double cos_theta = newGlycosidicProj.dot(oldGlycosidicProj);
        double sin_theta = n1.dot(newGlycosidicProj.cross(oldGlycosidicProj));
        double angle2 = atan2(sin_theta, cos_theta);
        Matrix3d R2 = AngleAxisd(angle2, n1).toRotationMatrix();

        // Combine rotation and translation
        Matrix3d R = R2 * R1;
        Vector3d translation = oldC1 - R * newC1;

        for (const auto& atom : baseTmpl) {
            if (atom.atomName == "C1'") {
                continue;
            }
            Atom newAtom = atom;
            Vector3d coord(atom.x, atom.y, atom.z);
            Vector3d newCoord = R * coord + translation;
            newAtom.atomType = "ATOM";
            newAtom.x = newCoord.x();
            newAtom.y = newCoord.y();
            newAtom.z = newCoord.z();
            newAtom.chainID = chainID[0];
            newAtom.resID = resID;
            newAtom.resName = newResName;
            newResidue.push_back(newAtom);
        }

        chain[resID] = newResidue;
        cout << "Residue " << chainID << resID << " has " << newResidue.size()
             << " atoms after replacement." << endl;
    }
}

void WritePDB(const Structure& structure, ostream& out) {
    int count = 0;
    for (const auto& chainPair : structure) {
        for (const auto& resPair : chainPair.second) {
            for (const auto& atom : resPair.second) {
                string element = atom.atomName.empty() ? " " : atom.atomName.substr(0, 1);
                out << setw(6) << left << atom.atomType
                    << right << setw(5) << ++count << " "
                    << left << setw(4) << atom.atomName << right << " "
                    << setw(3) << atom.resName
                    << " " << atom.chainID
                    << setw(4) << atom.resID << "    "
                    << fixed << setprecision(3)
                    << setw(8) << atom.x
                    << setw(8) << atom.y
                    << setw(8) << atom.z
                    << "  1.00  0.00           "
                    << setw(1) << element << "  \n";
            }
        }
        out << "TER\n";
    }
    out << "END\n";
}

bool HasHETATM(const Structure& structure) {
    const set<string> standardBases = {
        "A", "U", "G", "C", "DA", "DT", "DG", "DC"
    };
    for (const auto& [chainID, chain] : structure) {
        bool isNucleicAcid = false;
        for (const auto& [resID, residue] : chain) {
            if (!residue.empty() && standardBases.count(residue[0].resName)) {
                isNucleicAcid = true;
                break;
            }
        }
        if (!isNucleicAcid) {
            continue;
        }
        for (const auto& [resID, residue] : chain) {
            if (!residue.empty() && residue[0].atomType == "HETATM") {
                return true;
            }
        }
    }
    return false;
}

char GetAvailableChainID(const Structure& structure) {
    set<char> used;
    for (const auto& [chainID, _] : structure) {
        used.insert(chainID);
    }
    for (char c = 'A'; c <= 'Z'; ++c) {
        if (!used.count(c)) return c;
    }
    for (char c = 'a'; c <= 'z'; ++c) {
        if (!used.count(c)) return c;
    }
    return 'X';
}

string CreatePLIPPDB(const Structure& structure, char newChain, const string& filename) {
    const set<string> standardBases = {
        "A", "U", "G", "C", "DA", "DT", "DG", "DC"
    };

    Structure modified = structure;

    for (const auto& [chainID, chain] : structure) {
        bool isNucleicAcid = false;
        for (const auto& [resID, residue] : chain) {
            if (!residue.empty() && standardBases.count(residue[0].resName)) {
                isNucleicAcid = true;
                break;
            }
        }
        if (!isNucleicAcid) {
            continue;
        }

        for (const auto& [resID, residue] : chain) {
            if (!residue.empty() && residue[0].atomType == "HETATM") {
                for (const auto& atom : residue) {
                    Atom newAtom = atom;
                    newAtom.chainID = newChain;
                    newAtom.resID = atom.resID;
                    modified[newChain][atom.resID].push_back(newAtom);
                }
            }
        }
    }

    ofstream out(filename);
    if (!out) {
        return "";
    }
    WritePDB(modified, out);
    return filename;
}

int main(int argc, char** argv) {
    if (argc < 6 || ((argc - 2) % 4 != 0)) {
        cerr << "Usage: " << argv[0] << " <pdb_file> (<chainID> <resID_begin> <resID_end> <sequence>) ...\n";
        cerr << "Example 1: " << argv[0] << " examples/1emh.pdb B 2 6 ACGAT\n";
        cerr << "Example 2: " << argv[0] << " template.pdb A 1 8 ACGTACGT B 10 12 CGT" << endl;
        return 1;
    }

    fs::path inputPath = argv[1];
    Structure structure = ParsePDB(inputPath.string());
    Structure originalStructure = structure;

    if (structure.empty()) {
        cerr << "Error: No structure data found in " << inputPath.filename() << "." << endl;
        return 1;
    }

    fs::path dataPath = ResolveTemplatePath(argv[0]);
    if (!LoadBaseResidue(dataPath)) {
        cerr << "Error: Failed to load base residue templates." << endl;
        return 1;
    }

    string suffix;
    for (int i = 2; i < argc; i += 4) {
        string chainID = argv[i];
        int resID_begin = stoi(argv[i + 1]);
        int resID_end = stoi(argv[i + 2]);
        string seq = argv[i + 3];

        if (chainID.size() != 1) {
            cerr << "Error: Invalid chain ID '" << chainID << "' (must be a single character)." << endl;
            return 1;
        }
        if (resID_begin > resID_end || (resID_end - resID_begin + 1) != static_cast<int>(seq.size())) {
            cerr << "Error: Sequence length does not match the residue range for chain " << chainID << "." << endl;
            return 1;
        }

        suffix += "_" + chainID + "_" + to_string(resID_begin) + "_" + to_string(resID_end) + "_" + seq;
        HelixFromTmpl(structure, chainID, resID_begin, resID_end, seq);
    }

    fs::path outputPath = inputPath.stem().string() + suffix + inputPath.extension().string();
    ofstream outfile(outputPath);
    if (!outfile) {
        cerr << "Error: Failed to create output file " << outputPath << "!" << endl;
        return 1;
    }
    WritePDB(structure, outfile);
    outfile.close();
    cout << "Mutated PDB file generated as " << outputPath.filename() << "." << endl;

    string outdir = outputPath.string() + "_outputFolder";
    bool origHasHet = HasHETATM(originalStructure);
    bool mutHasHet = HasHETATM(structure);
    char copiedChain = (origHasHet || mutHasHet) ? GetAvailableChainID(structure) : '\0';

    bool plipOrigOk = false;
    if (origHasHet) {
        string tempPDB = CreatePLIPPDB(originalStructure, copiedChain, inputPath.string() + "_temp.pdb");
        if (tempPDB.empty()) {
            cerr << "Error creating temporary PDB file for PLIP." << endl;
            return 1;
        }
        plipOrigOk = runPlip(tempPDB, outdir + "/original");
        fs::remove(tempPDB);
    } else {
        plipOrigOk = runPlip(inputPath.string(), outdir + "/original");
    }

    bool plipMutOk = false;
    if (mutHasHet) {
        string tempPDB = CreatePLIPPDB(structure, copiedChain, outputPath.string() + "_temp.pdb");
        if (tempPDB.empty()) {
            cerr << "Error creating temporary PDB file for PLIP." << endl;
            return 1;
        }
        plipMutOk = runPlip(tempPDB, outdir + "/mutated");
        fs::remove(tempPDB);
    } else {
        plipMutOk = runPlip(outputPath.string(), outdir + "/mutated");
    }

    if (!plipOrigOk || !plipMutOk) {
        cerr << "Skipping interaction comparison because PLIP did not complete successfully." << endl;
        return 0;
    }

    InteractionMap oriMap = parseReport(outdir + "/original/report.xml");
    InteractionMap mutMap = parseReport(outdir + "/mutated/report.xml");

    compareInteractionsByType(oriMap, mutMap, copiedChain);
    writeInteractionsToMarkdown(oriMap, mutMap, copiedChain, outdir + "/report.md");

    return 0;
}
