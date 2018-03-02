import sys
import os
import yaml

vcf_tmpl = """\
##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=MENDACITY_MODES,Number=.,Type=String,Description="inheritance modes expected by mendacity">
##INFO=<ID=MENDACITY_NOT_MODES,Number=.,Type=String,Description="inheritance modes excluded by mendacity">
##contig=<ID=chr1,length=249250621,assembly=hg19>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s
"""

class Sample(object):
    __slots__ = ('id', 'paternal_id', 'maternal_id', 'sex', 'affected', 'alt_count')
    def __init__(self, values):
        self.id = values['id']
        self.affected = dict(affected='2', unaffected='1', carrier='1', unknown='-9')[values['status']]
        self.sex = dict(female='2', male='1', unknown='-9')[values['sex']]
        self.paternal_id = '-9' if values.get('father') is None else values['father']['id']
        self.maternal_id = '-9' if values.get('mother') is None else values['mother']['id']

    def __str__(self):
        return "fam_id\t%s\t%s\t%s\t%s\t%s" % (self.id, self.paternal_id, self.maternal_id, str(self.sex), self.affected)

    def __repr__(self):
        return "Sample(%s)" % self.id

    def tex(self, alt_count=1):
        return "\pstPerson[condition={aff} sex={sex} insidetext={gt}]{{{id}}}".format(aff="affected, " if self.affected == '2' else "normal", 
                sex={"2":"female", "1":"male", "-9":"unknown"}[self.sex],
                gt=["0/0", "0/1", "1/1", "./."][alt_count],
                id=self.id)

def gt(alt_count):
    ad = "%d,%d" % (10*(2-alt_count), 10*alt_count)
    fields = []
    fields.append(["0/0", "0/1", "1/1", "./."][alt_count])
    if alt_count in (0, 1, 2):
        fields.extend([ad, "20", "99", ["0,10,20", "10,0,20", "10,20,0"][alt_count]])

    return ":".join(fields)


class Pedigree(object):
    def __init__(self, fam_list):
        self.d = fam_list
        self.samples = [Sample(s) for s in fam_list]
        self.spouses = [(self[s.paternal_id], self[s.maternal_id]) for s in self.samples if not "-9" in (s.paternal_id, s.maternal_id)]

    def __getitem__(self, sample_id):
        for s in self.samples:
            if s.id == sample_id: return s
        raise KeyError(sample_id)

    @property
    def sample_ids(self):
        return [x.id for x in self.samples]


    def variant(self, alts, inheritance_modes, not_inheritance_modes, chrom="chr1", position="2345", ID="."):
        assert len(alts) == len(self.samples)
        sample_gts = [gt(a) for a in alts]

        info = []
        if inheritance_modes:
            info.append("MENDACITY_MODES=%s" % ",".join(inheritance_modes))
        if not_inheritance_modes:
            info.append("MENDACITY_NOT_MODES=%s" % ",".join(not_inheritance_modes))

        return "\t".join(map(str, ["chr1", position, ID, "T", "G", 50, "PASS", ";".join(info), "GT:AD:DP:GQ:PL"] + sample_gts))

    @property
    def tex(self):
        tmpl = """
\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{pstricks}
\usepackage{pst-pdgr}

\begin{document}
%s

\rput(0.5,1.5){\pstPerson[male, belowtext=1-1]{A}}
\rput(2.5,1.5){\pstPerson[affected, female, belowtext=1-2]{B}}
\rput(1.5,0.6){\pstPerson[male,belowtext=2-1]{C}}
\pstRelationship[descentnode=AB, rellinecmd=ncangle, angleA=90, angleB=90, descentnodepos=1.5, broken, brokenpos=1.2]{A}{B}
\ncline{AB}{C}


\end{document}
"""
    

    def pdf(self, prefix):
        rm(prefix + ".aux")
        rm(prefix + ".log")
        rm(prefix + ".pdf")
        with open(prefix + ".tex", "w") as fh:
            for s in self.samples:
                print s.tex()


def rm(path):
    try:
        os.unlink(path)
    except OSError:
        pass

def main(args=sys.argv[1:]):
    import os
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--prefix", default=None, help="prefix for output, default is basename of yaml")
    p.add_argument("yaml", help="path to YAML file containing tests and pedigree")

    a = p.parse_args(args)
    if a.prefix is None or a.prefix == ".":
        a.prefix = os.path.basename(a.yaml)
        if a.prefix.endswith(".yaml"):
            a.prefix = a.prefix[:-5]
        if a.prefix.endswith(".yml"):
            a.prefix = a.prefix[:-4]

    a.yaml = os.path.abspath(os.path.expandvars(os.path.expanduser(a.yaml)))

    y = yaml.load(open(a.yaml))
    p = Pedigree(y.pop('pedigree'))
    p.pdf(a.prefix)

    with open(a.prefix + ".ped", 'w') as fh:
        fh.write("\n".join(str(s) for s in p.samples))
        fh.write("\n")

        sys.stderr.write("wrote ped file to %s\n" % fh.name)

    cases = y.pop("cases")

    fh = open(a.prefix + ".vcf", 'w')
    fh.write(vcf_tmpl % "\t".join(p.sample_ids))
    position=55516888   # coding variant on chr1
    for i, case in enumerate(cases):
        for j, alt in enumerate(case['alts']):
            ID = "id%d_%d" % (i, j)
            fh.write(p.variant(alt, case.get('modes', []), case.get('not-modes', []), position=position, ID=ID))
            position += 4
            fh.write("\n")

    sys.stderr.write("wrote VCF file to %s\n" % fh.name)

if __name__ == "__main__":
    main()
