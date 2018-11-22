import os
import math
import pandas as pd
from citrination_client import CitrinationClient
from pypif import pif
from pypif.obj import *
from pif_updater.pore_statistics import *


def parse_csv(csv_file_dir):

    """
    Takes in csv file from dataset 73, returns pif system
    :return:
    """
    for f in os.listdir(csv_file_dir):
        if ".csv" in f:
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
                          'mean pore spacing', 'max pore diameter']

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

                system.properties.append(prop)

            print(pif.dumps(system.ids))
            outfile_path = csv_file_dir+f.replace('.csv', '.json')
            pif.dump(system, open(outfile_path, 'w'))


def get_files_from_dataset(dataset_id, download_path):

    client = CitrinationClient(os.environ['CITRINATION_ADAPT_API_KEY'], site='https://adapt.citrination.com')
    print(client)
    files = client.data.get_dataset_files(dataset_id=dataset_id)
    print(files)
    client.data.download_files(files, destination=download_path)


def breakout_pifs(base_input_dir):

    input_files =['P001_B001-nohough.json', 'P002_B001-nohough.json', 'P003_B001-nohough.json', 'P004_B001-nohough.json',
                  'P005_B001-nohough.json', 'P005_B002-nohough.json', 'P006_B001-nohough.json', 'P014_B001-nohough.json']

    for f in os.listdir(base_input_dir):

        if f in input_files:

            input_dir = base_input_dir+f.strip(".json")+"/"

            if not os.path.isdir(input_dir):
                os.makedirs(input_dir)

            if ".json" in f:
                systems = pif.load(open(base_input_dir+f))

                for system in systems:
                    for prep in system.preparation:
                        if prep.name == 'printing':
                            for det in prep.details:
                                if det.name == 'row':
                                    row_id = det.scalars
                                if det.name == 'column':
                                    column_id = det.scalars
                    if row_id >= 10:
                        sample_id = input_dir.split("/")[-2].split("-")[0]+"_"+column_id+str(row_id)
                    else:
                        sample_id = input_dir.split("/")[-2].split("-")[0]+"_"+column_id+"0"+str(row_id)

                    system.ids = [Id(name='Sample ID', value=sample_id)]

                    if "P001_B001" in f:
                        system.preparation.append(ProcessStep(name="Plate heat treatment", details=Value(name="Heat treatment performed", scalars="YES")))
                    else:
                        system.preparation.append(ProcessStep(name="Plate heat treatment", details=Value(name="Heat treatment performed", scalars="NO")))



                    outfile_path = input_dir+sample_id+".json"
                    pif.dump(system, open(outfile_path, "w"))
                    print("DUMPED: ", outfile_path)


def add_porosity_data_to_pif(csv_dir, pif_dir):

    for f in os.listdir(csv_dir):
        if ".json" in f:
            porosity_data_system = pif.load(open(csv_dir+f, "r"))

            sample_id = porosity_data_system.ids[0].value
            x = "_".join(sample_id.split('_')[0:2])

            main_system_path = pif_dir+x+"-nohough/"+sample_id+".json"

            try:
                main_system = pif.load(open(main_system_path, "r"))

                print(main_system_path)

                main_system.properties = porosity_data_system.properties

                pif.dump(main_system, open(main_system_path, 'w'))

            except Exception as e:
                print(e)


def merge_pifs(base_input_dir):

    for f in os.listdir(base_input_dir):

        if os.path.isdir(base_input_dir+f):

            print(f)
            systems = []

            for record in os.listdir(base_input_dir+f):
                if ".json" in record:
                    infile_path = base_input_dir+f+"/"+record
                    system = pif.load(open(infile_path, 'r'))
                    systems.append(system)

            outfile_path = base_input_dir+f+"_w-porosity.json"
            pif.dump(systems, open(outfile_path, 'w'))
            print('DUMPED: ', outfile_path, len(systems))


def add_classifier_to_pifs(base_input_dir):

    systems_with_props = []

    for f in os.listdir(base_input_dir):

        new_systems = []
        selected_prop_names = ['max pore diameter', 'mean pore diameter', 'fraction porosity',
                               'median pore spacing', 'median pore diameter']

        if "_w-porosity.json" in f:
            print(f)
            infile_path = base_input_dir+f
            old_systems = pif.load(open(infile_path, 'r'))
            print(len(old_systems))
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
                        if prop.name == 'max pore diameter':
                            mpd = float(prop.scalars.value)
                            if mpd > 200:
                                new_system.properties.append(Property(name='Pore size warning', scalars='RED'))
                            #elif 75 < mpd < 100:
                            #    new_system.properties.append(Property(name='Pore size warning', scalars='YELLOW'))
                            else:
                                new_system.properties.append(Property(name='Pore size warning', scalars='GREEN'))
                            new_system.properties.append(Property(name="log max pore diameter", scalars=Scalar(
                                value=math.log(max_pore_diameter(mpd)))))

                    systems_with_props.append(new_system)
                new_systems.append(new_system)
            outfile_path = infile_path.replace(".json", "_refined.json")
            pif.dump(new_systems, open(outfile_path, 'w'))

    outfile_path = base_input_dir+"all_pifs_with_props.json"
    pif.dump(systems_with_props, open(outfile_path, 'w'))


def upload_pifs(base_input_dir, dataset_id):

    for f in os.listdir(base_input_dir):
        if '-nohough.json' in f:
            print("UPLOADING: ", f)
            system = pif.load(open(base_input_dir+f, 'r'))
            print(base_input_dir+f)
            result = client.data.upload(dataset_id, base_input_dir+f, dest_path=f)
            print(result.__dict__)

if __name__ == "__main__":

    base_download_path = "/Users/cborg/Box Sync/Mines Open Lead [MOL]/projects/IN718/porosity/"
    client = CitrinationClient(os.environ['CITRINATION_ADAPT_API_KEY'], 'https://adapt.citrination.com')

    # get input csv files
    # get_files_from_dataset(dataset_id='74', download_path=base_download_path+"74/")

    # get pif files that need updating
    # get_files_from_dataset(dataset_id="73", download_path=base_download_path+"73/")

    # breakout records into individual pif files
    # breakout_pifs(base_input_dir=base_download_path+"73/")

    # parse csv
    # parse_csv(csv_file_dir=base_download_path+"74/")

    # Add data to corresponding pif file
    # add_porosity_data_to_pif(csv_dir=base_download_path+"74/", pif_dir=base_download_path+"73/")

    # re-merge pifs
    # merge_pifs(base_input_dir=base_download_path+"73/")

    # create 3 category prop for model
    add_classifier_to_pifs(base_input_dir=base_download_path+"73/")

    # upload to new dataset
    # upload_pifs(base_download_path+"73/", 78)