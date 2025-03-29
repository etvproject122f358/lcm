import webbrowser

def create_simbad_link(star_name):
    star_name = star_name.replace("+", "%2B")
    if " " in star_name:
        star_name = star_name.replace(" ", "+")
    base_url = f"https://simbad.cds.unistra.fr/simbad/sim-basic?Ident={star_name}&submit=SIMBAD+search"
    return base_url

def open_simbad_link(star_name):
    if star_name:
        simbad_link = create_simbad_link(star_name)
        webbrowser.open(simbad_link)
    else:
        print("Please enter a star name.")
