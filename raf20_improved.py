import requests
from shapely.geometry import Point
from shapely.ops import transform
import pyproj
import math

# ---------------------------------------------------------------------------
# Chargement de la grille RAF20
# ---------------------------------------------------------------------------

def charger_raf20(chemin='RAF20.tac'):
    """Charge la grille RAF20 depuis un fichier .tac et retourne
    le tableau 2D ainsi que les paramètres de la grille."""
    with open(chemin) as f:
        line = f.readline()
        params = line.split()

        min_lon  = float(params[0])
        max_lon  = float(params[1])
        min_lat  = float(params[2])
        max_lat  = float(params[3])
        step_lon = float(params[4])
        step_lat = float(params[5])

        # Chargement des données dans un tableau 2D
        # 0,0 = longitude minimale, latitude minimale
        # longitude croissante sur l'axe i, latitude croissante sur l'axe j
        h = int((max_lat - min_lat) / step_lat) + 1
        w = int((max_lon - min_lon) / step_lon) + 1
        raf20 = [[0 for _ in range(h)] for _ in range(w)]

        i = 0
        j = h - 1
        while line:
            line = f.readline()
            if not line.strip():
                continue
            h_elips_row = line.split()[::2]
            for h_elips in h_elips_row:
                c = [min_lon + i * step_lon, min_lat + j * step_lat, float(h_elips)]
                raf20[i][j] = c
                i += 1
                if i == w:
                    j -= 1
                    i = 0

    return raf20, min_lon, min_lat, step_lon, step_lat


def _interpoler(x_lambert93, y_lambert93, raf20, min_lon, min_lat, step_lon, step_lat):
    """Retourne la valeur T (ondulation du géoïde) interpolée en bilinéaire
    pour un point exprimé en Lambert 93 (EPSG:2154)."""
    lamb93 = pyproj.CRS('EPSG:2154')
    wgs84  = pyproj.CRS('EPSG:4326')
    to4326 = pyproj.Transformer.from_crs(lamb93, wgs84, always_xy=True).transform

    p_2154  = Point(x_lambert93, y_lambert93)
    p_4326  = transform(to4326, p_2154)

    n  = math.floor((p_4326.x - min_lon) / step_lon)
    q  = math.floor((p_4326.y - min_lat) / step_lat)
    fx = ((p_4326.x - min_lon) / step_lon) % 1
    fy = ((p_4326.y - min_lat) / step_lat) % 1

    T = (
        (1 - fx) * (1 - fy) * raf20[n][q][2]
        + fx      * (1 - fy) * raf20[n + 1][q][2]
        + (1 - fx) * fy      * raf20[n][q + 1][2]
        + fx       * fy      * raf20[n + 1][q + 1][2]
    )
    return T


# ---------------------------------------------------------------------------
# Fonctions publiques de conversion
# ---------------------------------------------------------------------------

def altitude_vers_ellipsoidale(x_lambert93, y_lambert93, z_altitude,
                                raf20, min_lon, min_lat, step_lon, step_lat):
    """Convertit une altitude terrain (NGF-IGN69) en hauteur ellipsoïdale.

    Paramètres
    ----------
    x_lambert93, y_lambert93 : coordonnées en Lambert 93 (EPSG:2154)
    z_altitude               : altitude terrain (m, NGF-IGN69)
    raf20, min_lon, …        : grille RAF20 chargée par charger_raf20()

    Retourne
    --------
    h_ellipsoidale : hauteur au-dessus de l'ellipsoïde WGS84 (m)
    T              : ondulation du géoïde (m)
    """
    T = _interpoler(x_lambert93, y_lambert93, raf20, min_lon, min_lat, step_lon, step_lat)
    h_ellipsoidale = z_altitude + T
    return h_ellipsoidale, T


def ellipsoidale_vers_altitude(x_lambert93, y_lambert93, h_ellipsoidale,
                                raf20, min_lon, min_lat, step_lon, step_lat):
    """Convertit une hauteur ellipsoïdale en altitude terrain (NGF-IGN69).

    Paramètres
    ----------
    x_lambert93, y_lambert93 : coordonnées en Lambert 93 (EPSG:2154)
    h_ellipsoidale           : hauteur au-dessus de l'ellipsoïde WGS84 (m)
    raf20, min_lon, …        : grille RAF20 chargée par charger_raf20()

    Retourne
    --------
    z_altitude : altitude terrain (m, NGF-IGN69)
    T          : ondulation du géoïde (m)
    """
    T = _interpoler(x_lambert93, y_lambert93, raf20, min_lon, min_lat, step_lon, step_lat)
    z_altitude = h_ellipsoidale - T
    return z_altitude, T

# ---------------------------------------------------------------------------
# Récupération des élévations d'un point (EPSG:2154) donné via l'API d'altimétrie de GeoPF
# https://cartes.gouv.fr/aide/fr/guides-utilisateur/utiliser-les-services-de-la-geoplateforme/calcul-altimetrique/
# ---------------------------------------------------------------------------
def get_elevation_from_lambert93(x_ref, y_ref):
    lamb93 = pyproj.CRS('EPSG:2154')
    wgs84 = pyproj.CRS('EPSG:4326')
    to4326 = pyproj.Transformer.from_crs(lamb93, wgs84, always_xy=True).transform
    p_2154 = Point(x_ref, y_ref)
    p_4326 = transform(to4326, p_2154)
    z_ref = get_elevation_from_wgs84(p_4326.x, p_4326.y)
    return z_ref

# ---------------------------------------------------------------------------
# Récupération des élévations d'un point (EPSG:4326) donné via l'API d'altimétrie de GeoPF
# https://cartes.gouv.fr/aide/fr/guides-utilisateur/utiliser-les-services-de-la-geoplateforme/calcul-altimetrique/
# ---------------------------------------------------------------------------
def get_elevation_from_wgs84(lon, lat):
    url = "https://data.geopf.fr/altimetrie/1.0/calcul/alti/rest/elevation.json"
    params = {
        "lon": lon,
        "lat": lat,
        "resource": "ign_rge_alti_wld"
    }
    response = requests.get(url, params=params)
    data = response.json()
    # Récupération de l'altitude
    z_ref = data["elevations"][0]["z"]
    return z_ref


# ---------------------------------------------------------------------------
# Test de round-trip
# ---------------------------------------------------------------------------

def test_round_trip(x_ref, y_ref, z_ref):
    """Vérifie que la conversion aller-retour retrouve bien les valeurs d'origine."""

    grille, min_lon, min_lat, step_lon, step_lat = charger_raf20()

    print(f"Coordonnées Lambert 93 : x = {x_ref} m, y = {y_ref} m")
    print(f"Altitude terrain (NGF-IGN69) : z = {z_ref:.4f} m\n")

    # Aller : altitude terrain → hauteur ellipsoïdale
    h_ellips, T = altitude_vers_ellipsoidale(
        x_ref, y_ref, z_ref, grille, min_lon, min_lat, step_lon, step_lat
    )
    print(f"Ondulation du géoïde T    = {T:.4f} m")
    print(f"Hauteur ellipsoïdale      = {h_ellips:.4f} m")
    print(f"Différence (h - z)        = {h_ellips - z_ref:.4f} m")

    # Retour : hauteur ellipsoïdale → altitude terrain
    z_retour, _ = ellipsoidale_vers_altitude(
        x_ref, y_ref, h_ellips, grille, min_lon, min_lat, step_lon, step_lat
    )
    print(f"Altitude terrain retrouvée = {z_retour:.4f} m (attendu : {z_ref} m)\n")

    # Vérification numérique
    assert abs(z_retour - z_ref) < 1e-9, (
        f"Échec du round-trip : {z_retour} ≠ {z_ref}"
    )
    print("✅ Test round-trip réussi : z retrouvé = z de départ")


if __name__ == '__main__':
    #x_ref = 730871   # x en Lambert 93
    #y_ref = 6336988  # y en Lambert 93
    #z_ref = 918.0    # altitude terrain (m)

    # Point sur la commune de Gibel
    x_ref = 590058  # x en Lambert 93
    y_ref = 6247486  # y en Lambert 93
    z_ref = get_elevation_from_lambert93(x_ref, y_ref)  # altitude terrain (m)

    test_round_trip(x_ref, y_ref, z_ref)
