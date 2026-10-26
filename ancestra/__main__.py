"""``python -m ancestra <command>``: ``run`` (core pipeline), ``eval`` (toy evaluation)."""
import sys


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: python -m ancestra {run,eval} ...\n"
              "  run   proteomes + guide tree -> HOGs (OrthoXML), see ancestra.pipeline\n"
              "  eval  compare predicted OrthoXML with simulated truth (ancestra.evaluate_toy)")
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "run":
        from .pipeline import main as m
    elif cmd == "eval":
        from .evaluate_toy import main as m
    else:
        print(f"unknown command {cmd!r}", file=sys.stderr)
        return 2
    return m(rest)


if __name__ == "__main__":
    sys.exit(main())
