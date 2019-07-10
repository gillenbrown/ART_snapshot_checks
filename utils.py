def print_and_write(info, file_obj, silent):
    if file_obj is not None:
        file_obj.write(str(info) + "\n")
    if not silent:
        print(info)