from multiprocessing import freeze_support
import traceback
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from sys import exit

def main():

    res = urlopen('http://just-the-time.appspot.com/')
    time_str = res.read().decode('utf-8').strip()

    date = time_str.split()[0]
    year, month, day = [int(s) for s in date.split("-")]

    print("\n*** ------------------------------------- ***")
    print("*** RGA: Reaction Generator and Analyzer  ***")
    print("*** Copyright @ Triet Le & Lam Huynh " + str(year) + " ***")
    print("*** ------------------------------------- ***\n")

    # Setting the expiry date: Nov 17, 2019

    expiry_day, expiry_month, expiry_year = 31, 12, 2020

    if year > expiry_year or (year == expiry_year and month > expiry_month) or \
            (year == expiry_year and month == expiry_month and day > expiry_day):
        print("Your version has been expired on December 31, 2020. Please contact the authors to get an updated version. Thank you !!!")
        input("Press Enter to continue...")
        exit(0)

    input("Press Enter to continue...")

if __name__ == '__main__':

    freeze_support()

    try:
        main()
    except (HTTPError, URLError) as e:

        print("\n*** ------------------------------------- ***")
        print("*** RGA: Reaction Generator and Analyzer  ***")
        print("*** Copyright @ Triet Le & Lam Huynh 2019 ***")
        print("*** ------------------------------------- ***\n")
        print("Please check your Internet connection or contact the authors for more information. Thank you !!!\n")

        with open("error.out", "w") as fout:
            fout.write("Please check your Internet connection or contact the authors for more information. Thank you !!!\n")

        input("Press Enter to continue...")
        exit(0)

    except Exception as e:

        print("There is something wrong with the program. Please check the details in error.out file. Thank you !!!\n")

        with open("error.out", "w") as fout:
            fout.write(str(traceback.format_exc()) + "\n")

        input("Press Enter to continue...")
        exit(0)
