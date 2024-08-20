import sys


class ProgressBar:
    def __init__(self, total, prefix="", suffix="", decimals=1, length=50, fill="█"):
        self.total = total
        self.prefix = prefix
        self.suffix = suffix
        self.decimals = decimals
        self.length = length
        self.fill = fill
        self.iteration = 0

    def update(self, iteration=None):
        self.iteration = iteration if iteration is not None else self.iteration + 1
        sys.stdout.write(str(self))
        sys.stdout.flush()

    def finish(self):
        self.update(self.total)
        print()

    def __str__(self):
        percent = f"{100 * (self.iteration / float(self.total)):.{self.decimals}f}"
        filled = int(self.length * self.iteration // self.total)
        bar = self.fill * filled + "-" * (self.length - filled)
        return f"\r{self.prefix} |{bar}| {percent}% {self.suffix}"

    @classmethod
    def create(cls, total, prefix="", suffix="", decimals=1, length=50, fill="█"):
        bar = cls(total, prefix, suffix, decimals, length, fill)
        return bar.update


class ProgresBarArray:
    """A class to manage multiple progress bars."""

    pass
    # TODO
