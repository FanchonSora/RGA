import time, os.path
from rdkit import Chem

speciesFile = input("Please enter the path of the species file: ")
if not os.path.exists(speciesFile):
    print("The path you input does not exist. Please try again !!!")
    input("Press Enter to continue...")
    exit(0)

found = False

speciesList = {}

start = time.perf_counter()

lineCnt = 1

print("Start checking the database of species...\n")

with open(speciesFile) as fin, open("SMILES_chk_res.out", "w") as fout:

    fout.write("The result of SMILES checking for the file " + speciesFile + ":\n")

    for line in fin:

        line = line.rstrip("\n")

        if line == "":
            continue

        # toks = line.split()
        toks = line.split(",")

        species = toks[0]

        mol = Chem.MolFromSmiles(species)

        lineCnt += 1

        if not mol:
            print(species, "is not a valid SMILES. Please check again !!!")
            print("------------------------")
            # fout.write(species + " is not a valid SMILES. Please check again !!!" + "\n")
            # fout.write("------------------------")

            fout.write(line + "\t")
            fout.write("!!Warning: invalid SMILES\n")

            found = True
            continue

        canonicalSmiles = Chem.MolToSmiles(mol, isomericSmiles=True)

        if canonicalSmiles != species:
            print("Please double-check the species", species)
            print("Suggested canonical SMILES fix:", canonicalSmiles)

            fout.write(line + "\t")
            fout.write("!!Warning: suggested SMILES: " + canonicalSmiles)

            if canonicalSmiles in speciesList:
                fout.write(" - Duplicate SMILES: " + canonicalSmiles + " on the line " + str(speciesList[canonicalSmiles]) + "\n")
                print(species, "is a duplicate of previous", canonicalSmiles, "on the line", str(speciesList[canonicalSmiles]))
            else:
                speciesList[canonicalSmiles] = lineCnt
                fout.write("\n")

            if not found:
                found = True

            print("------------------------\n")

            # fout.write("Please double-check the species " + species + "\n")
            # fout.write("Suggested canonical SMILES fix: " + canonicalSmiles + "\n")
            # fout.write("------------------------\n")
        else:
            fout.write(line)

            if canonicalSmiles in speciesList:
                print("Please double-check the species", species)
                print(species, "is a duplicate of previous", canonicalSmiles, "on the line", str(speciesList[canonicalSmiles]))
                print("------------------------\n")

                fout.write("\t!!Warning: Duplicate SMILES: " + canonicalSmiles + " on the line " + str(speciesList[canonicalSmiles]) + "\n")
            else:
                speciesList[canonicalSmiles] = lineCnt
                fout.write("\n")

if not found:
    print("All species are in correct form.")
    # fout.write("All species are in correct form.\n")

print("End checking the database of species...")

print("Please check the SMILES_chk_res.out for the results\n")

print("Execution time:", (time.perf_counter() - start), "s.")
print("--------------------------------\n")

input("Press Enter to continue...")
