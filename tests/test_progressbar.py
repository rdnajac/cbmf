from ..utils.progressbar import ProgressBar
import time


def test_progressbar():
    total = 100
    bar = ProgressBar(total, prefix="Progress:", suffix="Complete", length=40)

    for i in range(total + 1):
        bar.update(i)
        time.sleep(0.01)  # 10 ms wait

    bar.finish()


# if __name__ == "__main__":
#     test_progressbar()


def test_two_progressbars():
    total = 100
    bar1 = ProgressBar(total, prefix="Progress 1:", suffix="Complete", length=40)
    bar2 = ProgressBar(total, prefix="Progress 2:", suffix="Complete", length=40)

    for i in range(total + 1):
        bar1.update(i)
        bar2.update(i)
        time.sleep(0.01)  # 10 ms wait

    bar1.finish()
    bar2.finish()
