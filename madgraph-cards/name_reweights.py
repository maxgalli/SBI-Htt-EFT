import sys

def main(argv):
    inpath = argv[1]
    outpath = argv[2]

    outlines = []
    with open(inpath, "r") as infile:
        inlines = list(infile.readlines())
        for iline, line in enumerate(inlines):
            if line.startswith("launch"):
                next_lines = inlines[iline+1:iline+5]
                
                new_line = "launch --rwgt_name="
                for pline in next_lines:
                    pline_split = pline.split()
                    pname = pline_split[1]
                    pval = pline_split[2].replace(".", "p")

                    new_line += "{}_{}_".format(pname, pval)

                outlines.append(new_line[:-1] + "\n")
            else:
                outlines.append(line)

    with open(outpath, "w") as outfile:
        outfile.writelines(outlines)

    return 0

if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))