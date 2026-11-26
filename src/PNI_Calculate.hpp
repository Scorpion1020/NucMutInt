#pragma once

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "tinyxml2.h"

struct Clash {
    std::string mutatedAtom;
    std::string neighborAtom;
    double distance;
    double vdwSum;
};

struct Interaction {
    std::string proteinRes;
    std::string ligandRes;
    std::string distance;
    std::string angle;

    // Deduplicate and diff by residue pair so minor coordinate noise does not create false diffs.
    bool operator<(const Interaction& other) const {
        return std::tie(proteinRes, ligandRes) < std::tie(other.proteinRes, other.ligandRes);
    }

    std::string toString() const {
        std::string s = proteinRes + " <-> " + ligandRes;
        if (!distance.empty() && distance != "-") {
            s += " [dist: " + distance;
            if (!angle.empty() && angle != "-") {
                s += ", angle: " + angle;
            }
            s += "]";
        }
        return s;
    }
};

using InteractionMap = std::map<std::string, std::set<Interaction>>;

inline const std::vector<std::string> interactionTypes = {
    "hydrogen_bonds",
    "pi_stacks",
    "pi_cation_interactions",
    "salt_bridges",
    "hydrophobic_interactions",
    "water_bridges"
};

inline bool involvesChain(const Interaction& inter, char chainID) {
    if (chainID == '\0') {
        return false;
    }
    const std::string token = std::string(":") + chainID + ":";
    return inter.proteinRes.find(token) != std::string::npos ||
           inter.ligandRes.find(token) != std::string::npos;
}

inline InteractionMap parseReport(const std::string& xmlPath) {
    InteractionMap result;
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(xmlPath.c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Failed to load " << xmlPath << std::endl;
        return result;
    }

    tinyxml2::XMLElement* report = doc.FirstChildElement("report");
    if (!report) {
        return result;
    }

    for (tinyxml2::XMLElement* site = report->FirstChildElement("bindingsite");
         site;
         site = site->NextSiblingElement("bindingsite")) {
        tinyxml2::XMLElement* interactionsElem = site->FirstChildElement("interactions");
        if (!interactionsElem) {
            continue;
        }

        for (const std::string& type : interactionTypes) {
            tinyxml2::XMLElement* container = interactionsElem->FirstChildElement(type.c_str());
            if (!container) {
                continue;
            }

            for (tinyxml2::XMLElement* item = container->FirstChildElement();
                 item;
                 item = item->NextSiblingElement()) {
                auto getText = [item](const char* tag) -> const char* {
                    tinyxml2::XMLElement* elem = item->FirstChildElement(tag);
                    return (elem && elem->GetText()) ? elem->GetText() : nullptr;
                };

                const char* restype = getText("restype");
                const char* reschain = getText("reschain");
                const char* resnr = getText("resnr");

                const char* ligtype = getText("restype_lig");
                const char* ligchain = getText("reschain_lig");
                const char* lignr = getText("resnr_lig");

                if (!restype || !reschain || !resnr || !ligtype || !ligchain || !lignr) {
                    continue;
                }

                std::string distStr = "-";
                std::string angleStr = "-";

                if (type == "hydrogen_bonds") {
                    const char* dHA = getText("dist_h-a");
                    const char* dDA = getText("dist_d-a");
                    const char* ang = getText("don_angle");
                    if (dHA && dDA) {
                        distStr = std::string(dHA) + " A (H-A) / " + dDA + " A (D-A)";
                    } else if (dDA) {
                        distStr = std::string(dDA) + " A (D-A)";
                    }
                    if (ang) {
                        angleStr = std::string(ang) + " deg";
                    }
                } else if (type == "water_bridges") {
                    const char* dAW = getText("dist_a-w");
                    const char* dDW = getText("dist_d-w");
                    const char* ang = getText("don_angle");
                    if (dAW && dDW) {
                        distStr = std::string(dAW) + " A (A-W) / " + dDW + " A (D-W)";
                    }
                    if (ang) {
                        angleStr = std::string(ang) + " deg";
                    }
                } else if (type == "pi_stacks") {
                    const char* cd = getText("centdist");
                    const char* ang = getText("angle");
                    const char* stackType = getText("type");
                    if (cd) {
                        distStr = std::string(cd) + " A";
                    }
                    if (ang) {
                        angleStr = std::string(ang) + " deg";
                        if (stackType) {
                            angleStr += std::string(" (") + stackType + ")";
                        }
                    }
                } else {
                    const char* d = getText("dist");
                    if (d) {
                        distStr = std::string(d) + " A";
                    }
                }

                Interaction inter{
                    std::string(restype) + ":" + reschain + ":" + resnr,
                    std::string(ligtype) + ":" + ligchain + ":" + lignr,
                    distStr,
                    angleStr
                };
                result[type].insert(inter);
            }
        }
    }

    return result;
}

inline void compareInteractionsByType(const InteractionMap& ori,
                                      const InteractionMap& mut,
                                      char copiedChainID) {
    static const std::set<Interaction> emptySet;
    for (const auto& type : interactionTypes) {
        const auto& oriSet = ori.count(type) ? ori.at(type) : emptySet;
        const auto& mutSet = mut.count(type) ? mut.at(type) : emptySet;

        std::cout << "=== " << type << " ===" << std::endl;

        std::cout << "  Lost:" << std::endl;
        for (const auto& inter : oriSet) {
            if (!mutSet.count(inter) && !involvesChain(inter, copiedChainID)) {
                std::cout << "    " << inter.toString() << std::endl;
            }
        }

        std::cout << "  Gained:" << std::endl;
        for (const auto& inter : mutSet) {
            if (!oriSet.count(inter) && !involvesChain(inter, copiedChainID)) {
                std::cout << "    " << inter.toString() << std::endl;
            }
        }

        std::cout << std::endl;
    }
}

inline bool runPlip(const std::string& pdb, const std::string& outdir) {
    std::error_code ec;
    std::filesystem::create_directories(outdir, ec);

#ifdef _WIN32
    const std::string nullRedirect = " > NUL 2>&1";
#else
    const std::string nullRedirect = " > /dev/null 2>&1";
#endif

    const std::string cmd = "plip -f \"" + pdb + "\" -x -o \"" + outdir + "\"" + nullRedirect;
    std::cout << "[PLIP] Running: " << cmd << std::endl;
    const int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "[PLIP] Warning: command exited with code " << ret
                  << ". Ensure PLIP is installed and available on PATH." << std::endl;
        return false;
    }
    return true;
}

inline void writeInteractionsToMarkdown(const InteractionMap& ori,
                                        const InteractionMap& mut,
                                        char copiedChainID,
                                        const std::string& markdownOutputPath,
                                        const std::vector<Clash>& clashes = {}) {
    std::ofstream markdownFile(markdownOutputPath);
    if (!markdownFile.is_open()) {
        std::cerr << "Failed to open Markdown file for writing: " << markdownOutputPath << std::endl;
        return;
    }

    static const std::set<Interaction> emptySet;
    markdownFile << "# Interaction Difference Report\n\n";

    markdownFile << "## Steric Clash Check\n\n";
    if (clashes.empty()) {
        markdownFile << "_No severe steric clashes detected on mutated bases._\n\n";
    } else {
        markdownFile << "| Mutated Atom | Neighbor Atom | Distance (A) | VdW Sum (A) |\n";
        markdownFile << "| :--- | :--- | :---: | :---: |\n";
        markdownFile << std::fixed << std::setprecision(2);
        for (const auto& c : clashes) {
            markdownFile << "| `" << c.mutatedAtom << "` | `" << c.neighborAtom
                         << "` | " << c.distance << " | " << c.vdwSum << " |\n";
        }
        markdownFile << "\n";
    }

    markdownFile << "## Non-Covalent Interaction Changes\n\n";
    for (const auto& type : interactionTypes) {
        const auto& oriSet = ori.count(type) ? ori.at(type) : emptySet;
        const auto& mutSet = mut.count(type) ? mut.at(type) : emptySet;

        std::vector<const Interaction*> lost;
        std::vector<const Interaction*> gained;

        for (const auto& inter : oriSet) {
            if (!mutSet.count(inter) && !involvesChain(inter, copiedChainID)) {
                lost.push_back(&inter);
            }
        }
        for (const auto& inter : mutSet) {
            if (!oriSet.count(inter) && !involvesChain(inter, copiedChainID)) {
                gained.push_back(&inter);
            }
        }

        markdownFile << "### " << type << "\n\n";
        if (lost.empty() && gained.empty()) {
            markdownFile << "_No changes._\n\n";
            continue;
        }

        markdownFile << "| Status | Protein Residue | Ligand / Nucleic Acid | Distance | Angle |\n";
        markdownFile << "| :--- | :--- | :--- | :--- | :--- |\n";
        for (const auto* inter : lost) {
            markdownFile << "| Lost | `" << inter->proteinRes << "` | `" << inter->ligandRes
                         << "` | " << inter->distance << " | " << inter->angle << " |\n";
        }
        for (const auto* inter : gained) {
            markdownFile << "| Gained | `" << inter->proteinRes << "` | `" << inter->ligandRes
                         << "` | " << inter->distance << " | " << inter->angle << " |\n";
        }
        markdownFile << "\n";
    }

    markdownFile.close();
    std::cout << "Report successfully written to " << markdownOutputPath << std::endl;
}
