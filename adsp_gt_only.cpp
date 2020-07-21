#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

struct Variant {
	// the CHROM, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT columns in VCF
	std::string chr, id, ref, alt, qual, filter, info, format;
	int pos, ac, an; // position, AC, and AN
	std::vector<std::string> gt;
	std::vector<std::string> dp;
	std::vector<std::string> gq;
	bool pass;
};

//Get the basic info the first to nineth columns in VCF
inline void GetBasicInfo (Variant & var, std::stringstream & ss) {
	ss >> var.chr;
	ss >> var.pos;
	ss >> var.id;
	ss >> var.ref;
	ss >> var.alt;
	ss >> var.qual;
	ss >> var.filter;
	ss >> var.info;
	ss >> var.format;
}

// Keep falgs and ABHet in INFO and discard others
// Return True if the variant is good.
bool KeepFlag (std::string & info) {
	std::stringstream ss(info);
	std::string token, tmp;
	bool all_flagged = true, all_het_fail = true;
	while (std::getline(ss, token, ';')) { // Divid by ';'
		if (token.find ("VFLAGS_") != std::string::npos ||
			token.find ("ABHet_") != std::string::npos) { // Only keep element with VFLAGS_ and ABHet_
			if (!tmp.empty()) tmp += (';' + token); // Not the first element, needs ';'
			else tmp = token; // the first element in INFO

			if (token.find ("VFLAGS_") != std::string::npos) {
				std::size_t pos = token.find('=');
				if (token.substr(pos + 1) == "0")
					all_flagged = false;
			}

			if (token.find ("ABHet_") != std::string::npos) {
				std::size_t pos = token.find('=');
				float het = std::stod(token.substr(pos + 1));
				if (het > 0.25 && het < 0.75)
					all_het_fail = false;
			}
		}
	}

	info = tmp; // Replace the old info
	//return !(all_flagged & all_het_fail);
	return !all_flagged;
}

//Extract GT value of each sample and push them to var.gt
//If all GTs are missing, then return true.
bool ExtractGT (Variant & var, std::stringstream & ss) {
	//var.ac = 0;
	//var.an = 0;
	std::string token;
	//bool all_missing = true;
	bool first = true;
	while (std::getline(ss, token, '\t')) { // Divid by tab
		if (first) { // The first one is the empty one.
			first = false;
			continue;
		}
		std::size_t pos = token.find(':'); // GT is the first info in ADSP VCF
		const std::string gt = token.substr(0, pos); // Extract GT
		var.gt.push_back(gt); // Store GT

		/*
		all_missing &= (gt == "./.");
	
		if (gt[0] != '.') var.an++;
		if (gt[2] != '.') var.an++;
		if (gt[0] != '.' && gt[0] != '0') var.ac++;
		if (gt[2] != '.' && gt[2] != '0') var.ac++;
		*/

		/*
		pos = token.find(':'), pos + 1; // AD
		std::size_t pos1 = token.find((':'), pos + 1); // DP
		var.dp.push_back(token.substr(pos + 1, pos1 - pos - 1));
		pos = pos1;
		pos1 = token.find((':'), pos + 1); // GQ
		var.gq.push_back(token.substr(pos + 1, pos1 - pos - 1));
		*/
		
	}
	//return all_missing;
	return false;
}

void PrintVariant(const Variant & var) {
	//if (!var.pass) return;
	std::cout << var.chr << "\t"
		<< var.pos << "\t"
		<< var.id << "\t"
		<< var.ref << "\t"
		<< var.alt << "\t"
		<< var.qual << "\t"
		<< var.filter << "\t"
		<< var.info << "\t"
		//<< "AC=" << var.ac << ";AF=" << var.ac/var.an << ";AN=" << var.an << ";" << var.info << "\t"
		<< var.format;
	// Print GT of each sample
	for (unsigned int i = 0; i < var.gt.size(); ++i) // The first one GT is empty, so we skip it.
		std::cout << "\t" << var.gt[i];
	std::cout << std::endl;
}

int main (int argc, char** argv) {

	// USAGE: zcat vcf.gz | <PROGRAM> | gzip > output.vcf.gz

	std::ios_base::sync_with_stdio(false); // Fastern IO
	std::cin.tie(NULL); // Fastern IO

	std::string line; // buffer
	while (std::getline(std::cin, line)) { // Catch a line from stdin
		if (line[0] == '#') { // for header section, write out directly
			std::cout << line << std::endl;
		} else {
			std::stringstream ss(line); // Convert line to stringstream
			Variant var;
			GetBasicInfo(var, ss);
			var.format = "GT"; // Only keep GT
			//var.pass = KeepFlag(var.info); // rephase info; keep only VFLAGS and ABHet
			ExtractGT(var, ss);
			PrintVariant(var);
		}
	}

	return 0;
}
