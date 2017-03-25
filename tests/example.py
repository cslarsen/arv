from arv import load, unphased_match as match

genome = load("genome.txt")

print("You are a {gender} with {color} eyes and {complexion} skin.".format(
    gender     = "man" if genome.y_chromosome else "woman",
    complexion = "light" if genome["rs1426654"] == "AA" else "dark",
    color      = match(genome["rs12913832"], {"AA": "brown",
                                              "AG": "brown or green",
                                              "GG": "blue"})))
