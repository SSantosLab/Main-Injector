import argparse
import sys





def hello_world(args=None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument('-f', type=str)

    args = p.parse_args(args)

    if args.f:
        print(args.f)