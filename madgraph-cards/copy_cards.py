import sys
import os
import shutil

def main(argv):
    assert len(argv) == 3, "Wrong number of arguments"

    inpath = argv[1]
    outpath = argv[2]
    outname = outpath.split("/")[-1] if outpath[-1] != "/" else outpath.split("/")[-2]
    inname = inpath.split("/")[-1] if inpath[-1] != "/" else inpath.split("/")[-2]

    os.mkdir(outpath)
    card_names = os.listdir(inpath)

    for card in card_names:
        out_card = card.replace(inname, outname)
        shutil.copy(inpath + "/" + card, outpath + "/" + out_card)

    return 0

if __name__ == "__main__":
    print("\nFinished with exit code:", main(sys.argv))