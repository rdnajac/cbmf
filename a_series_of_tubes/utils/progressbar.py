import sys

def print_progress_bar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=50,
    fillchar="█",

):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled = int(length * iteration // total)
    progress = fillchar * filled + "-" * (length - filled)

    progress_bar = f'\r{prefix} |{progress}| {percent}% {suffix}'
    
    # Print the progress bar and end with carriage return
    sys.stdout.write(progress_bar)
    sys.stdout.flush()

    # Print a new line once the progress bar is complete
    if iteration == total:
        # print a complete prog 
        progress=100
        filled = length
        progress = fillchar * length 
        progress_bar = f'\r{prefix} |{progress}| {percent}% {suffix}\n'
        sys.stdout.write(progress_bar)
        sys.stdout.flush()


# def create_progress_bar(
#     total, prefix="", suffix="", decimals=1, length=50, fill="█", print_end="\r"
# ):
#     def print_progress_bar(iteration):
#         percent = ("{0:." + str(decimals) + "f}").format(
#             100 * (iteration / float(total))
#         )
#         filled_length = int(length * iteration // total)
#         bar = fill * filled_length + "-" * (length - filled_length)
#         pr.print_color(
#             f"\r{prefix} |{bar}| {percent}% {suffix}", fg=Color.CYAN, end=print_end
#         )

#     return print_progress_bar


# class ProgressBar:
#     def __init__(self, total, prefix="", suffix="", decimals=1, length=50, fill="█"):
#         self.progress_bar = create_progress_bar(
#             total, prefix, suffix, decimals, length, fill
#         )

#     def update(self, iteration):
#         self.progress_bar(iteration)

#     def finish(self):
#         print()  # New line after progress bar is complete
