#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

struct Variant {
	// the CHROM, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT columns in VCF
	std::string chr, id, ref, alt, qual, filter, info, format;
	int pos, ac, an; // position, AC, and AN.
	int qc_dp, qc_gq, qc_both; // DP < 10 and/or GQ <20.
	int ac_qc, an_qc; // AC and AN after QC.
	std::vector<std::string> gt;
	std::vector<std::string> dp;
	std::vector<std::string> gq;
	bool pass;

	bool clean() {
		chr.clear();
		id.clear();
		ref.clear();
		alt.clear();
		qual.clear();
		filter.clear();
		info.clear();
		format.clear();
		pos=0, ac=0, an=0;
		qc_dp=0, qc_gq=0, qc_both=0;
		ac_qc=0, an_qc=0;
		gt.clear();
		dp.clear();
		gq.clear();
		pass = false;
	}
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

//Add "ori" as a prefix to AN, AC, AF, and DP in INFO field.
inline void AddPrefix (const std::string & prefix, const std::string & target, std::string & info) {
	std::size_t found = info.find(target);
	if (found != std::string::npos)
		info.insert(found, prefix);
	else
		std::cerr << "Warning: Cannot find " << target << " in " << info << std::endl;
		
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

//Extract GT,DP and GQ values of each sample and push them to var.gt, var.dp, and var.gq.
bool ExtractInfo (Variant & var, std::stringstream & ss) {
	//var.ac = 0;
	//var.an = 0;
	std::string token;
	//bool all_missing = true;
	bool first = true;
	const std::string missingGt = "./.";
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
		*/

		int an = 0, ac = 0;

		if (gt[0] != '.') an++;
		if (gt[2] != '.') an++;
		if (gt[0] != '.' && gt[0] != '0') ac++;
		if (gt[2] != '.' && gt[2] != '0') ac++;

		var.an += an;
		var.ac += ac;
		var.an_qc += an;
		var.ac_qc += ac;
		
		// For AD
		pos = token.find(':'), pos + 1; // AD

		// For DP
		std::size_t pos1 = token.find((':'), pos + 1); // DP
		//var.dp.push_back(token.substr(pos + 1, pos1 - pos - 1));
		const bool qc_dp_fail = (std::atoi(token.substr(pos + 1, pos1 - pos - 1).c_str()) < 10) ? true : false;
		if (qc_dp_fail) {
			var.qc_dp++;
			var.an_qc -= an;
			var.ac_qc -= ac;
			*(var.gt.rbegin()) = missingGt;
		}

		// For GQ
		pos = pos1;
		pos1 = token.find((':'), pos + 1); // GQ
		//var.gq.push_back(token.substr(pos + 1, pos1 - pos - 1));
		const bool qc_gq_fail = (std::atoi(token.substr(pos + 1, pos1 - pos - 1).c_str()) < 20) ? true : false;
		if (qc_gq_fail) var.qc_gq++;
		if (qc_gq_fail && !qc_dp_fail) {
			var.an_qc -= an;
			var.ac_qc -= ac;
			*(var.gt.rbegin()) = missingGt;
		}

		if (qc_dp_fail && qc_gq_fail) var.qc_both++;
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
		//<< var.info << "\t"
		<< "AC=" << var.ac_qc << ";AF=" << var.ac_qc/var.an_qc << ";AN=" << var.an_qc << ";" << var.info << "\t"
		<< var.format;
	// Print GT of each sample
	for (unsigned int i = 0; i < var.gt.size(); ++i) // The first one GT is empty, so we skip it.
		std::cout << "\t" << var.gt[i];
	std::cout << std::endl;
}

int main (int argc, char** argv) {

	// USAGE: zcat vcf.gz | <PROGRAM> | gzip > output.vcf.gz 2> info.txt

	std::ios_base::sync_with_stdio(false); // Fastern IO
	std::cin.tie(NULL); // Fastern IO

	std::string line; // buffer
	std::cerr << "AC\tAN\tQC_AC\tQC_AN\tDP<10\tGQ<20\tDP_GQ" << std::endl;
	while (std::getline(std::cin, line)) { // Catch a line from stdin
		if (line[0] == '#') { // for header section, write out directly
			std::cout << line << std::endl;
		} else {
			std::stringstream ss(line); // Convert line to stringstream
			Variant var;
			var.clean();
			GetBasicInfo(var, ss);
			AddPrefix("ori_", "AC=", var.info);
			AddPrefix("ori_", "AN=", var.info);
			AddPrefix("ori_", "AF=", var.info);
			AddPrefix("ori_", "DP=", var.info);
			var.format = "GT"; // Only keep GT
			//var.pass = KeepFlag(var.info); // rephase info; keep only VFLAGS and ABHet
			ExtractInfo(var, ss); // Change low-qual genotypes to missing
			PrintVariant(var); // Output vcf to stdout
			std::cerr << var.ac << "\t" << var.an << "\t" << var.ac_qc << "\t" 
			<< var.an_qc << "\t" << var.qc_dp << "\t" << var.qc_gq << "\t" << var.qc_both << std::endl;
		}
	}

	return 0;
}
