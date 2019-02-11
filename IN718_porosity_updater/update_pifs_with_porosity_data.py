import os
import math
import sys
import pandas as pd
from citrination_client import CitrinationClient
from pypif import pif
from pypif.obj import *
from IN718_porosity_updater.pore_statistics import *
sys.path.insert(0, '/Users/cborg/projects/community_projects/')
from community_projects.pycc_utils import pycc_wrappers


def parse_csv(csv_file_dir, pif_dir):

    """
    Takes in csv file from dataset 73, returns pif system
    _full.csv = total volume of part
    :return:
    """
    csv_files = [f for f in os.listdir(csv_file_dir) if ".csv" in f and "_full" not in f]
    full_csv_files = [f for f in os.listdir(csv_file_dir) if "_full.csv" in f]

    for f in csv_files:
        df = pd.read_csv(csv_file_dir+f, encoding="utf-16")

        system = ChemicalSystem()
        sample_id = f.strip(".csv")
        system.ids = [Id(name='Sample ID', value=sample_id)]

        cm_x = [Scalar(value=x) for x in df['Center Of Mass X (µm)']]
        cm_y = [Scalar(value=x) for x in df['Center Of Mass Y (µm)']]
        cm_z = [Scalar(value=x) for x in df['Center Of Mass Z (µm)']]

        method = Method(name='porosity', software=Software(name='tracr', version='beta'))
        prop_x = Property(name='center of mass X', scalars=cm_x, units='$\mu m$', method=method)
        prop_y = Property(name='center of mass Y', scalars=cm_y, units='$\mu m$', method=method)
        prop_z = Property(name='center of mass Z', scalars=cm_z, units='$\mu m$', method=method)
        system.properties = [prop_x, prop_y, prop_z]

        # calc pore stats
        pore_stats = ['neighbor pore distance', 'median pore diameter', 'median pore spacing',
                      'mean pore spacing', 'max pore diameter',  'pore volume', 'pore diameters',
                      'stdev of pore diameters', 'total pores']

        for prop_name in pore_stats:
            prop = Property()
            prop.name = prop_name
            if prop_name == 'median pore diameter':
                prop.scalars = Scalar(value=median_pore_diameter(df['Volume (µm³)']))
                prop.units = "$\mu m$"
            if prop_name == 'neighbor pore distance':
                prop.scalars = [Scalar(value=x) for x in nearest_neighbor_distance(df['Center Of Mass X (µm)'],
                                                                                   df['Center Of Mass Y (µm)'],
                                                                                   df['Center Of Mass Z (µm)'])]
                prop.units = '$\mu m$'
            if prop_name == 'median pore spacing':
                prop.scalars = Scalar(value=median_pore_spacing(df['Center Of Mass X (µm)'],
                                                                df['Center Of Mass Y (µm)'],
                                                                df['Center Of Mass Z (µm)']))
                prop.units = '$\mu m$'
            if prop_name == 'mean pore spacing':
                prop.scalars = Scalar(value=mean_pore_spacing(df['Center Of Mass X (µm)'],
                                                              df['Center Of Mass Y (µm)'],
                                                              df['Center Of Mass Z (µm)']))
                prop.units = '$\mu m$'
            if prop_name == 'max pore diameter':
                prop.scalars = Scalar(value=max_pore_diameter(df['Volume (µm³)']))
                prop.units = '$\mu m$'

            if prop_name == 'pore volume':
                prop.scalars = [Scalar(value=x) for x in df['Volume (µm³)']]
                prop.units = '${\mu m}^3$'

            if prop_name == 'pore diameters':
                prop.scalars = [Scalar(value=x) for x in sphere_equivalent_diameter(df['Volume (µm³)'])]
                prop.units = '$\mu m$'

            if prop_name == 'stdev of pore diameters':
                prop.scalars = Scalar(value=round(np.std(sphere_equivalent_diameter(df['Volume (µm³)'])), 3))
                prop.units = '$\mu m$'

            if prop_name == 'total pores':
                prop.scalars = Scalar(value=len(df['Volume (µm³)']))

            system.properties.append(prop)

        print(pif.dumps(system.ids))
        outfile_path = pif_dir+f.replace('.csv', '.json')
        pif.dump(system, open(outfile_path, 'w'))

    for f in full_csv_files:

        df = pd.read_csv(csv_file_dir+f, encoding="utf-16")

        outfile_path = f.replace('_full.csv', '.json')

        if outfile_path in os.listdir(pif_dir):
            system = pif.load(open(pif_dir+outfile_path, 'r'))
            # system.properties.append(Property(name='Full part volume', scalars=df['Volume (µm³)'], units='${\mu m}^3$'))
            for prop in system.properties:
                if prop.name == 'pore volume':
                    total_porosity_vol = sum([sca.value for sca in prop.scalars])
                    fractional_porosity = round(float(total_porosity_vol / df['Volume (µm³)']), 6)
                    system.properties.append(Property(name='fraction porosity', scalars=fractional_porosity))

            pif.dump(system, open(pif_dir+outfile_path, 'w'))
            print("Fraction porosity calc: ", outfile_path)


def get_files_from_dataset(dataset_id, download_path):

    client = CitrinationClient(os.environ['CITRINATION_ADAPT_API_KEY'], site='https://adapt.citrination.com')
    print(client)
    files = client.data.get_dataset_files(dataset_id=dataset_id)
    print("Downloading {} files...".format(len(files)))
    client.data.download_files(files, destination=download_path)


def add_identifiers_to_pifs(systems, f):

    for system in systems:
        for prep in system.preparation:
            if prep.name == 'printing':
                for det in prep.details:
                    if det.name == 'row':
                        row_id = det.scalars
                    if det.name == 'column':
                        column_id = det.scalars
        if row_id >= 10:
            sample_id = f.replace("-nohough.json", "") + "_" + column_id + str(row_id)
        else:
            sample_id = f.replace("-nohough.json", "") + "_" + column_id + "0" + str(row_id)

        system.ids = [Id(name='Sample ID', value=sample_id)]

    return systems


def add_heat_treatment_to_pifs(systems, f):

    for system in systems:
        if "P001_B001" in f:
            system.preparation.append(
                ProcessStep(name="Plate heat treatment", details=Value(name="Heat treatment performed", scalars="YES")))
        else:
            system.preparation.append(
                ProcessStep(name="Plate heat treatment", details=Value(name="Heat treatment performed", scalars="NO")))
    return systems


def modify_master_dataset(master_branch_dir, develop_branch_dir):

    for f in os.listdir(master_branch_dir):

        if ".json" in f:
            systems = pif.load(open(master_branch_dir + f))
            systems = add_identifiers_to_pifs(systems, f)
            systems = add_heat_treatment_to_pifs(systems, f)
            systems = add_porosity_data_to_pifs(systems, base_download_path+"data/porosity_jsons/")
            systems = add_porosity_stats_to_pifs(systems)
            systems = add_pore_diameter_bucket_prop(systems)
            systems = remove_unverified_pore_data(systems)
            outfile_path = develop_branch_dir+f
            pif.dump(systems, open(outfile_path, "w"))
            print("DUMPED: ", outfile_path)

def remove_unverified_pore_data(systems):

    unverified_ids = ['P001_B001_X13', 'P001_B001_B03', 'P001_B001_B14']
    if system in systems:
        if system.ids[0].value in unverified_ids:
            system.properties = []

    return systems

def add_pore_diameter_bucket_prop(systems):

    for system in systems:
        if system.properties:
            for prop in system.properties:
                if prop.name == 'pore diameters':
                    system.properties.append(Property(name='pore diameter < 50 um', scalars=len([i for i in prop.scalars if float(i.value) < 50])))
                    system.properties.append(Property(name='pore diameter 50 < x < 100 um', scalars=len([i for i in prop.scalars if 50 < float(i.value) < 100])))
                    system.properties.append(Property(name='pore diameter 100 < x < 150 um', scalars=len([i for i in prop.scalars if 100 < float(i.value) < 150])))
                    system.properties.append(Property(name='pore diameter 150 < x < 200 um', scalars=len([i for i in prop.scalars if 150 < float(i.value) < 200])))
                    system.properties.append(Property(name='pore diameter x > 200 um', scalars=len([i for i in prop.scalars if float(i.value) > 200])))
    return systems

def add_porosity_data_to_pifs(systems, data_porosity_jsons):

    for f in os.listdir(data_porosity_jsons):
        if ".json" in f:
            porosity_data_system = pif.load(open(data_porosity_jsons+f, "r"))
            porosity_system_sample_id = porosity_data_system.ids[0].value
            for system in systems:
                main_system_sample_id = system.ids[0].value
                if main_system_sample_id == porosity_system_sample_id:
                    system.properties = porosity_data_system.properties

    return systems


def add_porosity_stats_to_pifs(systems):

    for system in systems:
        if system.properties:
            prop_names = [prop.name for prop in system.properties]
            for prop in system.properties:
                if prop.name == 'pore volume':

                    pore_volumes = [float(sca.value) for sca in prop.scalars]
                    r_squared_norm = round(qq_normal(sphere_equivalent_diameter(pore_volumes))[2]**2, 4)
                    r_squared_lognorm = round(qq_lognormal(sphere_equivalent_diameter(pore_volumes))[2]**2, 4)
                    system.properties.append(Property(name='r_squared_norm', scalars=r_squared_norm))
                    system.properties.append(Property(name='r_squared_lognorm', scalars=r_squared_lognorm))
                    if r_squared_norm > r_squared_lognorm:
                        system.properties.append(Property(name='dist_best_fit', scalars='NORM'))
                    else:
                        system.properties.append(Property(name='dist_best_fit', scalars='LOGNORM'))

                    if 'pore diameters' not in prop_names:
                        pore_diameters = [Scalar(value=x) for x in sphere_equivalent_diameter(pore_volumes)]
                        system.properties.append(Property(name='pore diameters', scalars=pore_diameters, units='$\mu m$'))

                    if 'stdev of pore diameters' not in prop_names:
                        stdev = Scalar(value=round(np.std(sphere_equivalent_diameter(pore_volumes)), 3))
                        system.properties.append(Property(name='stdev of pore diameters', scalars=stdev, units='$\mu m$'))

                    if 'total pores' not in prop_names:
                        total_pores = Scalar(value=len(pore_volumes))
                        system.properties.append(Property(name='total pores', scalars=total_pores))

                if prop.name == 'max pore diameter':
                    mpd = float(prop.scalars.value)
                    if mpd > 200:
                        system.properties.append(Property(name='Pore size warning', scalars='RED'))
                    else:
                        system.properties.append(Property(name='Pore size warning', scalars='GREEN'))

                    if mpd > 200:
                        system.properties.append(Property(name='Pore size warning (ternary)', scalars='RED'))
                    elif 75 < mpd < 200:
                        system.properties.append(Property(name='Pore size warning (ternary)', scalars='YELLOW'))
                    else:
                        system.properties.append(Property(name='Pore size warning (ternary)', scalars='GREEN'))

                    system.properties.append(Property(name="log max pore diameter", scalars=Scalar(value=math.log10(mpd))))

                if prop.name == 'median pore diameter':
                    if prop.scalars.value > 22:
                        system.properties.append(Property(name='Median pore classifier', scalars='>22 um'))
                    else:
                        system.properties.append(Property(name='Median pore classifier', scalars='<22 um'))

    return systems


def refine_to_relevant_props(develop_branch_dir, feature_branch_dir):

    porosity_props = ['max pore diameter', 'mean pore diameter', 'fraction porosity', 'median pore spacing',
                           'median pore diameter', 'log max pore diameter', 'Pore size warning',
                           'Pore size warning (ternary)', 'total pores', 'stdev of pore diameters', 'dist_best_fit',
                           'r_squared_norm', 'r_squared_lognorm', 'pore diameter < 50 um', 'pore diameter 50 < x < 100 um',
                           'pore diameter 100 < x < 150 um', 'pore diameter 150 < x < 200 um', 'pore diameter x > 200 um',
                           'Median pore classifier']

    mechanical_props = ['elastic modulus', 'elastic onset', 'yield strength', 'yield strain', 'ultimate strength',
                        'necking onset', 'fracture strength', 'total elongation', 'ductility', 'toughness']

    selected_prop_names = porosity_props + mechanical_props

    for f in os.listdir(develop_branch_dir):

        if ".json" in f:

            new_systems = []

            infile_path = develop_branch_dir + f
            old_systems = pif.load(open(infile_path, 'r'))
            print(infile_path, len(old_systems))
            for old_system in old_systems:
                new_system = ChemicalSystem()
                new_system.names = old_system.names
                new_system.references = old_system.references
                new_system.ids = old_system.ids
                new_system.preparation = old_system.preparation
                # new_system.sub_systems = old_system.sub_systems
                new_system.properties = []
                if old_system.properties:
                    for prop in old_system.properties:
                        if prop.name in selected_prop_names:
                            new_system.properties.append(prop)

                # mechanical props stored in subsystem
                if old_system.sub_systems:
                    for sub_system in old_system.sub_systems:
                        if sub_system.properties:
                            for prop in sub_system.properties:
                                if prop.name in selected_prop_names:
                                    new_system.properties.append(prop)

                new_systems.append(new_system)

            outfile_path = feature_branch_dir+f.replace(".json", "_refined.json")
            pif.dump(new_systems, open(outfile_path, 'w'))



def upload_pifs(base_input_dir, dataset_id):

    for f in os.listdir(base_input_dir):
        if ".json" in f:
            print("UPLOADING: ", base_input_dir+f)
            result = client.data.upload(dataset_id, base_input_dir+f, dest_path=f)
            print(result.__dict__)


def remove_outliers(base_input_dir):

    for f in os.listdir(base_input_dir):

        if "_refined.json" in f:
            print(f)
            infile_path = base_input_dir+f
            systems = pif.load(open(infile_path, 'r'))

            for system in systems:
                for prop in system.properties:
                    if prop.name == 'max pore diameter':
                        mpd = float(prop.scalars.value)
                        if mpd > 120:
                            print(pif.dumps(prop))
                            prop.scalars = ""

            outfile_path = infile_path.replace(".json", "_no_outliers.json")
            pif.dump(systems, open(outfile_path, 'w'))


# refines unlabeled records in design space to just records from P005_B002
def refine_design_space(input_dir, output_dir):

    for f in os.listdir(input_dir):

        if ".json" in f:
            infile_path = input_dir + f
            systems = pif.load(open(infile_path, 'r'))
            print(f, len(systems))

            refined_systems = []

            for system in systems:
                if not system.properties and "P005_B002" not in f:
                    pass
                else:
                    refined_systems.append(system)

            print(f, len(refined_systems))
            outfile_path = output_dir + f
            pif.dump(refined_systems, open(outfile_path, 'w'))


def refine_by_id(input_dir, output_dir):

    ids = ["P005_B002_V09", "P005_B002_U09", "P005_B002_W09", "P005_B002_O04", "P005_B002_P04", "P005_B002_V07",
           "P005_B002_Y09", "P005_B002_V04", "P005_B002_T09", "P005_B002_V10", "P005_B002_V06", "P005_B002_V08",
           "P005_B002_V02", "P005_B002_V01", "P005_B002_L06", "P005_B002_V05", "P005_B002_O03", "P005_B002_L07",
           "P005_B002_X10", "P005_B002_C14", "P005_B002_V11", "P005_B002_B14", "P005_B002_A15", "P005_B002_O02"]

    for f in os.listdir(input_dir):

        if ".json" in f:
            infile_path = input_dir + f
            systems = pif.load(open(infile_path, 'r'))
            print(f, len(systems))

            refined_systems = []

            for system in systems:
                if system.ids[0].value in ids:
                    system.properties = []
                refined_systems.append(system)

            print(f, len(refined_systems))
            outfile_path = output_dir + f
            pif.dump(refined_systems, open(outfile_path, 'w'))


if __name__ == "__main__":

    base_download_path = "/Users/cborg/Box Sync/Mines Open Lead [MOL]/projects/NAVSEA/IN718/"
    client = CitrinationClient(os.environ['CITRINATION_ADAPT_API_KEY'], 'https://adapt.citrination.com')

    # get input csv files from a data branch
    # get_files_from_dataset(dataset_id='74', download_path=base_download_path+"data/porosity_csvs/")

    # since there is no outfacing ingester, ingest csvs files
    # parse_csv(csv_file_dir=base_download_path+"data/porosity_csvs/", pif_dir=base_download_path+"data/porosity_jsons/")

    # get pif files from master branch
    # get_files_from_dataset(dataset_id="73", download_path=base_download_path+"master/LPBF_Inconel_718/")

    # modify master branch, pushes modified dataset to a develop branch
    # modify_master_dataset(master_branch_dir=base_download_path+"master/LPBF_Inconel_718/", develop_branch_dir=base_download_path+"develop/LPBF_Inconel_718/")

    refine_to_relevant_props(develop_branch_dir=base_download_path+"develop/LPBF_Inconel_718/", feature_branch_dir=base_download_path+"feature/IN718_refined_with_mech_props/")

    # upload to new dataset
    # upload_pifs(base_download_path+"develop/LPBF_Inconel_718/", 78)

    # remove_outliers(base_input_dir=base_download_path+"73/")
    # pif_check(base_download_path+"feature/IN718_refined/", base_download_path+"feature/ml_ready/")

    # refine_design_space(input_dir=base_download_path+"feature/IN718_refined/", output_dir=base_download_path+"feature/IN718_refined_design/")
    # refine_by_id(input_dir=base_download_path+"feature/IN718_refined_design/", output_dir=base_download_path+"feature/IN718_refined_design_week1/")