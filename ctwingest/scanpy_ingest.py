
from ctwingest.io import make_dir_or_complain, write_all_worksheet, delete_dir
from ctwingest.marker_table import run_pipe
import ctwingest.scanpyapi as ad_obj
import tarfile


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def scanpy_ingest(ad, worksheet_name, cluster_name, celltype_key):
    print("Attempt to gather cell type mapping")
    mapping = ad_obj.celltype_mapping(ad, cluster_name, celltype_key)

    if mapping is None:
        print("No Cell Type Mapping Found")
    else:
        print("Mapping Found: preview:", mapping.head())

    use_raw = ad_obj.has_raw(ad)
    xys = ad_obj.get_xys(ad, key="X_umap")

    print("Running marker generation")

    clustering = ad_obj.get_obs(ad, cluster_name)


    markers_df = run_pipe(ad, cluster_name)


    exp = ad_obj.get_expression(ad, use_raw)

    # Make the directory to tar up later.
    make_dir_or_complain(worksheet_name)
    write_all_worksheet(worksheet_name, xys=xys, exp=exp,
        clustering=clustering, markers=markers_df, celltype=mapping)

    ctw_filename = "%s.ctw.tgz" % worksheet_name
    make_tarfile(ctw_filename, source_dir=worksheet_name)
    delete_dir(worksheet_name)
