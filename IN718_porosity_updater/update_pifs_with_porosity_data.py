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
    :return:
    """
    for f in os.listdir(csv_file_dir):
        if ".csv" in f and "_full" not in f:
            df = pd.read_csv(csv_file_dir+f, encoding="utf-16")

            cm_x = [Scalar(value=x) for x in df['Center Of Mass X (µm)']]
            cm_y = [Scalar(value=x) for x in df['Center Of Mass Y (µm)']]
            cm_z = [Scalar(value=x) for x in df['Center Of Mass Z (µm)']]

            system = ChemicalSystem()
            sample_id = f.strip(".csv")
            system.ids = [Id(name='Sample ID', value=sample_id)]

            method = Method(name='porosity', software=Software(name='tracr', version='beta'))
            prop_x = Property(name='center of mass X', scalars=cm_x, units='$\mu m$', method=method)
            prop_y = Property(name='center of mass Y', scalars=cm_y, units='$\mu m$', method=method)
            prop_z = Property(name='center of mass Z', scalars=cm_z, units='$\mu m$', method=method)
            system.properties = [prop_x, prop_y, prop_z]

            # calc pore stats
            pore_stats = ['neighbor pore distance', 'median pore diameter', 'median pore spacing',
                          'mean pore spacing', 'max pore diameter',  'pore volume',
                          'stdev of pore diameter distribution', 'total pores']

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

                if prop_name == 'stdev of pore diameter distribution':
                    prop.scalars = Scalar(value=round(np.std(sphere_equivalent_diameter(df['Volume (µm³)'])), 3))
                    prop.units = '$\mu m$'

                if prop_name == 'total pores':
                    prop.scalars = Scalar(value=len(df['Volume (µm³)']))

                system.properties.append(prop)

            print(pif.dumps(system.ids))
            outfile_path = pif_dir+f.replace('.csv', '.json')
            pif.dump(system, open(outfile_path, 'w'))


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
            outfile_path = develop_branch_dir+f
            pif.dump(systems, open(outfile_path, "w"))
            print("DUMPED: ", outfile_path)


def add_porosity_data_to_pifs(systems, data_porosity_jsons):

    for f in os.listdir(data_porosity_jsons):

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
            for prop in system.properties:
                if prop.name == 'max pore diameter':
                    mpd = float(prop.scalars.value)
                    if mpd > 200:
                        system.properties.append(Property(name='Pore size warning', scalars='RED'))
                    else:
                        system.properties.append(Property(name='Pore size warning', scalars='GREEN'))

                    if mpd > 200:
                        system.properties.append(Property(name='Pore size warning (ternary)', scalars='RED'))
                    elif 75 < mpd < 100:
                        system.properties.append(Property(name='Pore size warning (ternary)', scalars='YELLOW'))
                    else:
                        system.properties.append(Property(name='Pore size warning (ternary)', scalars='GREEN'))

                    system.properties.append(Property(name="log max pore diameter", scalars=Scalar(value=math.log(max_pore_diameter(mpd)))))
    return systems


def refine_to_relevant_props(develop_branch_dir, feature_branch_dir):

    new_systems = []
    selected_prop_names = ['max pore diameter', 'mean pore diameter', 'fraction porosity', 'median pore spacing',
                           'median pore diameter', 'log max pore diameter', 'Pore size warning',
                           'Pore size warning (ternary)']

    for f in os.listdir(develop_branch_dir):

        if ".json" in f:

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

    # refine_to_relevant_props(develop_branch_dir=base_download_path+"develop/LPBF_Inconel_718/", feature_branch_dir=base_download_path+"feature/IN718_refined/")

    # Add data to corresponding pif file
    # add_porosity_data_to_pif(csv_dir=base_download_path+"74/", pif_dir=base_download_path+"73/")


    # upload to new dataset
    # upload_pifs(base_download_path+"develop/LPBF_Inconel_718/", 78)

    # remove_outliers(base_input_dir=base_download_path+"73/")