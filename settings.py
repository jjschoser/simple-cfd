from file_handler import *


def write_settings(settings_fname, init_header_fname, final_time, 
                   bc_lo=None, bc_hi=None, gamma=None, cfl=None, 
                   sdf_header_fname=None, out_header_base_fname=None, 
                   out_interval=None, vdb_base_fname=None, vdb_interval=None, 
                   vdb_start_idx=None):
    lines = ["init_header_fname = " + init_header_fname + "\n",
             "final_time = " + str(final_time) + "\n"]
    if bc_lo is not None:
        lines.append("bc_lo = " + " ".join([str(x) for x in bc_lo]) + "\n")
    if bc_hi is not None:
        lines.append("bc_hi = " + " ".join([str(x) for x in bc_hi]) + "\n")
    if gamma is not None:
        lines.append("gamma = " + str(gamma) + "\n")
    if cfl is not None:
        lines.append("cfl = " + str(cfl) + "\n")
    if sdf_header_fname is not None:
        lines.append("sdf_header_fname = " + sdf_header_fname + "\n")
    if out_header_base_fname is not None:
        lines.append("out_header_base_fname = " + out_header_base_fname + "\n")
    if out_interval is not None:
        lines.append("out_interval = " + str(out_interval) + "\n")
    if vdb_base_fname is not None:
        lines.append("vdb_base_fname = " + vdb_base_fname + "\n")
    if vdb_interval is not None:
        lines.append("vdb_interval = " + str(vdb_interval) + "\n")
    if vdb_start_idx is not None:
        lines.append("vdb_start_idx = " + str(vdb_start_idx) + "\n")
    create_dir_for_file(settings_fname)
    with open(settings_fname, "w") as f:
        f.writelines(lines)
