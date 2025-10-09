import os
import string
import random
import textwrap

from rich import print

def print_motd() -> None:
    """
    Prints a motd message with a randomized greeting and a unique identifier
    """

    from ..constants import APP_VERSION
    from time import gmtime, strftime

    width: int = 60

    print("[b]" + width * '~' + "[/b]")
    print("[b]|" + f"Biomolecule Annotation Toolkit v{APP_VERSION}".center(width - 2, " ") + "|[/b]")
    print("[b]|" + "".center(width - 2, " ") + "|[/b]")
    print("[b]|" + f" working directory: {os.getcwd()}".ljust(width - 2, " ") + "|[/b]")
    curr_time: str = strftime("%d-%m-%Y %H:%M:%S", gmtime())
    print("[b]|" + f" running at: {curr_time}".ljust(width - 2, " ") + "|[/b]")
    print("[b]" + width * '~' + "[/b]")
    print()

    for line in textwrap.wrap(
        f'"{random_quote().strip()}"', width, break_long_words=True
    ):
        print("[i]" + line + "[/i]")
    print()


def random_quote(from_file: str = "src/quotes.txt") -> str:
    """
    If quote file exists, returns one random quote from it. Otherwise, returns default one.
    :param from_file: path to the file where one line contains one quote
    :return: one random quote
    """
    if not os.path.exists(from_file):
        return "It is impossible to begin to learn that which one thinks one already knows. - Epictetus"

    all_quotes: list[str] = []
    with open(from_file, "r") as fd:
        for line in fd:
            if not line.startswith("#"):
                all_quotes.append(line)
    return random.choice(all_quotes)

def random_string(length: int) -> str:
    """
    Generates a random string of letters and numbers
    :param length: of the string
    :return: random string
    """
    return ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits)
                   for _ in range(length))


from .. import core_config
uid_length: int = core_config.getint("miscellaneous", "uid_length")


def random_uid() -> str:
    """
    :return: random unique identifier string with fixed length from the config
    """
    return random_string(uid_length)