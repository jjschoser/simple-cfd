from pathlib import Path

# Contents of this file were written with the help of ChatGPT

def get_data_fname(header_fname):
    path = Path(header_fname)
    return str(path.with_name(path.stem + "Data.dat"))


def remove_path(fname):
    path = Path(fname)
    return path.name


def add_path(src_fname, dst_fname):
    src_path = Path(src_fname)
    new_path = src_path.parent / dst_fname
    return str(new_path)


def remove_step_counter(fname):
    path = Path(fname)
    stem = path.stem
    i = len(stem) - 1
    while i >= 0 and stem[i].isdigit():
        i -= 1
    new_stem = stem[:i]
    return str(path.with_stem(new_stem))


def create_dir_for_file(fname):
    path = Path(fname)
    path.parent.mkdir(parents=True, exist_ok=True)


def get_header_fname_list(base_name, dir="."):
    header_files = []
    for path in Path(dir).glob(f"{base_name}*"):
        suffix = path.stem[len(base_name)+1:]
        if suffix.isdigit():
            num = int(suffix)
            header_files.append((num, path))
    header_files.sort(key=lambda x: x[0])
    return [file[1] for file in header_files]


def get_last_header_fname(base_name, dir="."):
    last_num = 0
    last_file = base_name + str(last_num)
    for path in Path(dir).glob(f"{base_name}*"):
        suffix = path.stem[len(base_name)+1:]
        if suffix.isdigit():
            num = int(suffix)
            if num > last_num:
                last_num = num
                last_file = path
    return last_file
