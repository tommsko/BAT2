import json
import os
from collections import defaultdict
from enum import Enum
from typing import Type, Tuple
from rich import print
from tabulate import tabulate


class StatsAggregator:
    def __init__(self, manager: 'StatsAggregatorManager') -> None:
        """
        Initializes (generic) stats aggregator
        :param manager: aggregator manager
        """
        self.manager: StatsAggregatorManager = manager

    def run(self):
        """
        Executes the stats aggregation and prints results to stdout.
        :return: None
        """
        raise NotImplementedError


class LoadStatsAggregator(StatsAggregator):
    """
    Aggregates statistics about how many simulations were loaded successfully
    """

    def __init__(self, manager: 'StatsAggregatorManager') -> None:
        super().__init__(manager)

    def run(self) -> None:
        sims_count: int = len(self.manager.simulation_folders)
        sims_loaded_count: int = self.manager._successfully_loaded

        load_type: dict[str, int] = defaultdict(int)
        for res in self.manager.result_files.values():
            if "simulation" in res and "load_strategy" in res["simulation"]:
                load_type[res["simulation"]["load_strategy"]] += 1
            else:
                load_type["failed"] += 1

        print("\n\n\n")
        print("[b][cyan]Simulation loading statistics [/cyan][/b]")
        print("[b]Shows how simulations were loaded into the workflow[/b]")
        print()

        print(f"Simulations provided: {sims_count}")
        print(f"[green]... of which successfully loaded: {sims_loaded_count}[/green]")
        print(f"[red]... of which unsuccessfully loaded: {sims_count - sims_loaded_count}[/red]")
        print()

        print("Counts of successfully loaded simulations vs load strategies:")
        table_content: list[list[str]] = [
            [type_str, count] for type_str, count in load_type.items()
        ]
        table_headers: list[str] = ["Load strategy", "Count"]
        print(tabulate(table_content, table_headers, tablefmt="simple_grid"))


class SignatureStatsAggregator(StatsAggregator):
    """
    Aggregates statistics about how molecules were described with signatures
    """

    def __init__(self, manager: 'StatsAggregatorManager') -> None:
        super().__init__(manager)

        self._all_signature_generators: list[str] = []
        self._analytical_signature_generators: list[str] = []
        self._descriptive_signature_generators: list[str] = []
        self._fallback_signature_generators: list[str] = []
        self._populate_signature_generators()

        self._all_molecules_cnt: int = 0
        self._sign_analytical: int = 0
        self._sign_descriptive_not_analytical: int = 0
        self._sign_fallback_not_descriptive: int = 0
        self._sign_none: int = 0

        self._signatures_attempted: dict[str, int] = defaultdict(int)
        self._signatures_generated: dict[str, int] = defaultdict(int)

        self._signatures_failed_to_generate: dict[str, set[str]] = defaultdict(set)

    def _populate_signature_generators(self) -> None:
        """
        Populates the lists of signature generators
        :return: None
        """
        from ..signatures.signature_generators.signature_generator import SignatureGeneratorType
        from ..signatures.signature_manager import SignatureType

        for signature_type in SignatureType:
            self._all_signature_generators.append(signature_type.name)
            match (signature_type.value)().generator_type:
                case SignatureGeneratorType.ANALYTICAL:
                    self._analytical_signature_generators.append(signature_type.name)
                case SignatureGeneratorType.DESCRIPTIVE:
                    self._descriptive_signature_generators.append(signature_type.name)
                case SignatureGeneratorType.FALLBACK:
                    self._fallback_signature_generators.append(signature_type.name)
                case _:
                    assert False, f"Unsupported signature type: {signature_type.name}"

    def _aggregate_stats(self) -> None:
        """
        Aggregates statistics about signatures from result files
        :return:
        """
        for res_name, res in self.manager.result_files.items():
            if "signatures" not in res:
                continue
            for segment_name in res["signatures"]:

                signatures_attempted: list[str] = res["signatures"][segment_name]["attempted"]
                signatures_generated: dict[str, str] = res["signatures"][segment_name]["generated"]

                self._all_molecules_cnt += 1

                if any(signature_type in self._analytical_signature_generators
                       for signature_type in signatures_generated.keys()):
                    self._sign_analytical += 1
                elif any(signature_type in self._descriptive_signature_generators
                       for signature_type in signatures_generated.keys()):
                    self._sign_descriptive_not_analytical += 1
                elif any(signature_type in self._fallback_signature_generators
                       for signature_type in signatures_generated.keys()):
                    self._sign_fallback_not_descriptive += 1

                for signature_type in signatures_attempted:
                    self._signatures_attempted[signature_type] += 1
                for signature_type in signatures_generated.keys():
                    self._signatures_generated[signature_type] += 1

                for signature_type in signatures_attempted:
                    if signature_type not in signatures_generated.keys():
                        self._signatures_failed_to_generate[signature_type].add(res_name)

        self._sign_none = self._all_molecules_cnt - (self._sign_analytical
                                                     + self._sign_descriptive_not_analytical
                                                     + self._sign_fallback_not_descriptive)

    def _get_signature_generator_table(self) -> Tuple[list[list[str]], list[str]]:
        """
        Generates attempted/generated signatures by signature generator table
        :return: table content and table headers
        """
        header: list[str] = ["Signature generator", "Input", "->", "Attempted at", "%", "Generated at", "%"]
        table_content: list[list[str]] = []
        for i, signature_type in enumerate(self._all_signature_generators):
            table_content.append([
                signature_type,
                str(self._all_molecules_cnt) if i == 0 else "",
                "->",
                str(self._signatures_attempted[signature_type]),
                str(round(self._signatures_attempted[signature_type]/self._all_molecules_cnt*100, 2)) + "%",
                str(self._signatures_generated[signature_type]),
                str(round(self._signatures_generated[signature_type] / self._all_molecules_cnt * 100, 2)) + "%",
            ])
        return table_content, header

    def run(self) -> None:
        self._aggregate_stats()

        print("\n\n\n")
        print("[b][cyan]Signature generation statistics [/cyan][/b]")
        print("[b]Shows how molecules were described into signatures[/b]")
        print()

        print(f"Molecules analyzed:                          {self._all_molecules_cnt} (100%)")
        print(f"[green]... w/ analytical signature:                 {self._sign_analytical} "
              f"({round(self._sign_analytical/self._all_molecules_cnt*100, 2)}%)[/green]")
        print(f"[green]... w/ descriptive signature w/o analytical: {self._sign_descriptive_not_analytical} "
              f"({round(self._sign_descriptive_not_analytical/self._all_molecules_cnt*100, 2)}%)[/green]")
        print(f"[green]... w/ fallback signature w/o descriptive:   {self._sign_fallback_not_descriptive} "
              f"({round(self._sign_fallback_not_descriptive/self._all_molecules_cnt*100, 2)}%)[/green]")
        print(f"[red]... w/o any signature:                       {self._sign_none} "
              f"({round(self._sign_none/self._all_molecules_cnt*100, 2)}%)[/red]")
        print()

        print("Performance of signature generators:")
        print(tabulate(*self._get_signature_generator_table(), tablefmt="simple_grid"))
        print()

        print("[magenta][bold]Export simulations where signature generator was attempted but failed?[/bold][/magenta]")
        if input("Enter 'y' to export, any other key to skip: ").lower() != 'y':
            return
        file_name: str = input("Enter the file name to export (default: signatures_failed_to_generate.json): ")
        if not file_name: file_name = "signatures_failed_to_generate.json"

        with open(file_name, 'w') as fd:
            json.dump({key: list(value) for key, value in self._signatures_failed_to_generate.items()}, fd, indent=4)


class IdentifierStatsAggregator(StatsAggregator):
    """
    Aggregates statistics about how signatures were identified
    """

    class IField(Enum):
        ATTEMPTED = 1
        IDENTITY = 2
        SIMILARITY = 3

    def __init__(self, manager: 'StatsAggregatorManager') -> None:
        super().__init__(manager)

        self._all_signature_generators: list[str] = []
        self._analytical_signature_generators: list[str] = []
        self._descriptive_signature_generators: list[str] = []
        self._fallback_signature_generators: list[str] = []
        self._populate_signature_generators()

        self._all_molecules_cnt: int = 0
        self._identified_identity: int = 0
        self._identified_similarity: int = 0
        self._identified_none: int = 0

        self._all_signatures_cnt: int = 0
        self._sign_analytical_ident: int = 0
        self._sign_descriptive_ident: int = 0
        self._sign_fallbacl_ident: int = 0
        self._sign_ident: int = 0
        self._sign_analytical_sim: int = 0
        self._sign_descriptive_sim: int = 0
        self._sign_fallbacl_sim: int = 0
        self._sign_sim: int = 0
        self._sign_unresolved: int = 0

        # signature_generator -> identifier_resolver -> field -> value
        self._identifier_stats: defaultdict[str, defaultdict[str, defaultdict[IdentifierStatsAggregator.IField, int]]] \
            = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        # signature_generator -> cnt
        self._signature_stats: defaultdict[str, int] = defaultdict(int)

        # results_file_name -> segment_name -> generated_signatures
        self._failed_to_identify: defaultdict[str, defaultdict[str, list[str]]] = defaultdict(lambda: defaultdict(list))

    def _populate_signature_generators(self) -> None:
        """
        Populates the lists of signature generators
        :return: None
        """
        from ..signatures.signature_generators.signature_generator import SignatureGeneratorType
        from ..signatures.signature_manager import SignatureType

        for signature_type in SignatureType:
            self._all_signature_generators.append(signature_type.name)
            match (signature_type.value)().generator_type:
                case SignatureGeneratorType.ANALYTICAL:
                    self._analytical_signature_generators.append(signature_type.name)
                case SignatureGeneratorType.DESCRIPTIVE:
                    self._descriptive_signature_generators.append(signature_type.name)
                case SignatureGeneratorType.FALLBACK:
                    self._fallback_signature_generators.append(signature_type.name)
                case _:
                    assert False, f"Unsupported signature type: {signature_type.name}"

    def _aggregate_stats(self) -> None:
        """
        Aggregates statistics about identifiers from result files
        :return:
        """
        for res_name, res in self.manager.result_files.items():
            if "signatures" not in res:
                continue

            segments: list[str] = res["signatures"].keys()
            for segment in segments:
                self._all_molecules_cnt += 1
                ident_identity: bool = False
                ident_similarity: bool = False

                if "identifiers" not in res or segment not in res["identifiers"]:
                    continue

                for signature_generator, identification_data in res["identifiers"][segment].items():
                    self._all_signatures_cnt += 1

                    self._signature_stats[signature_generator] += 1
                    for identifier_resolver in identification_data["attempted"]:
                        self._identifier_stats[signature_generator][identifier_resolver][IdentifierStatsAggregator.IField.ATTEMPTED] += 1

                    if identification_data["identity"] is not None:
                        ident_identity = True
                        self._identifier_stats[signature_generator][identification_data["identifier_name"]][IdentifierStatsAggregator.IField.IDENTITY] += 1
                        if signature_generator in self._analytical_signature_generators:
                            self._sign_analytical_ident += 1
                        elif signature_generator in self._descriptive_signature_generators:
                            self._sign_descriptive_ident += 1
                        elif signature_generator in self._fallback_signature_generators:
                            self._sign_fallbacl_ident += 1
                        else:
                            raise RuntimeError(f"Unknown signature generator type for '{signature_generator}'")

                    elif identification_data["similarity"] is not None:
                        ident_similarity = True
                        self._identifier_stats[signature_generator][identification_data["identifier_name"]][IdentifierStatsAggregator.IField.SIMILARITY] += 1
                        if signature_generator in self._analytical_signature_generators:
                            self._sign_analytical_sim += 1
                        elif signature_generator in self._descriptive_signature_generators:
                            self._sign_descriptive_sim += 1
                        elif signature_generator in self._fallback_signature_generators:
                            self._sign_fallbacl_sim += 1
                        else:
                            raise RuntimeError(f"Unknown signature generator type for '{signature_generator}'")

                    else:
                        self._sign_unresolved += 1

                if ident_identity:
                    self._identified_identity += 1
                elif ident_similarity:
                    self._identified_similarity += 1
                else:
                    self._identified_none += 1
                    self._failed_to_identify[res_name][segment] = list(res["signatures"][segment]["generated"].keys())

        self._sign_ident = self._sign_analytical_ident + self._sign_descriptive_ident + self._sign_fallbacl_ident
        self._sign_sim = self._sign_analytical_sim + self._sign_descriptive_sim + self._sign_fallbacl_sim


    def _get_identifier_resolution_table(self) -> Tuple[list[list[str]], list[str]]:
        """
        Generates attempted/generated identifiers by resolver table
        :return: table content and table headers
        """
        header: list[str] = ["Signature/identifier", "Input", "Attempted at", "%", "Primary identity at", "%", "Primary similarity at", "%"]
        table_content: list[list[str]] = []
        for signature_type, data in self._identifier_stats.items():
            signatures_cnt: int = self._signature_stats[signature_type]

            sum_ident: int = 0
            sum_sim: int = 0
            rows = []
            for identifier_name, stats in data.items():
                attempted_at: int = stats[IdentifierStatsAggregator.IField.ATTEMPTED]
                attempted_at_p: str = str(round(attempted_at / signatures_cnt * 100, 2))
                primary_identity_at: int = stats[IdentifierStatsAggregator.IField.IDENTITY]
                sum_ident += primary_identity_at
                primary_identity_at_p: str = str(round(primary_identity_at / signatures_cnt * 100, 2))
                primary_similarity_at: int = stats[IdentifierStatsAggregator.IField.SIMILARITY]
                sum_sim += primary_similarity_at
                primary_similarity_at_p: str = str(round(primary_similarity_at / signatures_cnt * 100, 2))
                rows.append([
                    f" ⤷ {identifier_name}",
                    "",
                    str(attempted_at),
                    attempted_at_p,
                    str(primary_identity_at),
                    primary_identity_at_p,
                    str(primary_similarity_at),
                    primary_similarity_at_p,
                ])

            sum_ident_p = str(round(sum_ident / signatures_cnt * 100, 2))
            sum_sim_p = str(round(sum_sim / signatures_cnt * 100, 2))

            table_content.append([
                signature_type,
                str(signatures_cnt),
                "",
                "",
                str(sum_ident),
                sum_ident_p,
                str(sum_sim),
                sum_sim_p,
            ])

            table_content.extend(rows)

            table_content.append(["", "", "", "", "", "", "", ""])

        return table_content, header

    def run(self) -> None:
        self._aggregate_stats()

        print("\n\n\n")
        print("[b][cyan]Identifier resolution statistics [/cyan][/b]")
        print("[b]Shows how molecules were resolved into identifiers[/b]")
        print()

        print(f"Molecules analyzed:                          {self._all_molecules_cnt} (100%)")
        print(f"[green] ⤷ resolved by exact-match:                  {self._identified_identity} ({round(self._identified_identity / self._all_molecules_cnt * 100, 2)}%)[/green]")
        print(f"[yellow] ⤷ resolved by similarity:                   {self._identified_similarity} ({round(self._identified_similarity / self._all_molecules_cnt * 100, 2)}%)[/yellow]")
        print(f"[red] ⤷ not resolved:                             {self._identified_none} ({round(self._identified_none / self._all_molecules_cnt * 100, 2)}%)[/red]")
        print()

        print(f"Signatures analyzed:                         {self._all_signatures_cnt} (100%)")
        print(f"[green] ⤷ exact-match resolution:                   {self._sign_ident} ({round(self._sign_ident / self._all_signatures_cnt * 100, 2)}%)[/green]")
        print(f"[green]    ⤷ analytical:                            {self._sign_analytical_ident} ({round(self._sign_analytical_ident / self._all_signatures_cnt * 100, 2)}%)[/green]")
        print(f"[green]    ⤷ descriptive:                           {self._sign_descriptive_ident} ({round(self._sign_descriptive_ident / self._all_signatures_cnt * 100, 2)}%)[/green]")
        print(f"[green]    ⤷ fallback:                              {self._sign_fallbacl_ident} ({round(self._sign_fallbacl_ident / self._all_signatures_cnt * 100, 2)}%)[/green]")
        print(f"[yellow] ⤷ similarity resolution:                    {self._sign_sim} ({round(self._sign_sim / self._all_signatures_cnt * 100, 2)}%)[/yellow]")
        print(f"[yellow]    ⤷ analytical:                            {self._sign_analytical_sim} ({round(self._sign_analytical_sim / self._all_signatures_cnt * 100, 2)}%)[/yellow]")
        print(f"[yellow]    ⤷ descriptive:                           {self._sign_descriptive_sim} ({round(self._sign_descriptive_sim / self._all_signatures_cnt * 100, 2)}%)[/yellow]")
        print(f"[yellow]    ⤷ fallback:                              {self._sign_fallbacl_sim} ({round(self._sign_fallbacl_sim / self._all_signatures_cnt * 100, 2)}%)[/yellow]")
        print(f"[red] ⤷ not resolved:                             {self._sign_unresolved} ({round(self._sign_unresolved / self._all_signatures_cnt * 100, 2)}%)[/red]")
        print()

        print("Performance of identifier resolvers")
        print(tabulate(*self._get_identifier_resolution_table(), tablefmt="simple_grid"))
        print()

        print("[magenta][bold]Export molecules without successful identification?[/bold][/magenta]")
        if input("Enter 'y' to export, any other key to skip: ").lower() != 'y':
            return
        file_name: str = input("Enter the file name to export (default: molecules_failed_to_identify.txt): ")
        if not file_name: file_name = "molecules_failed_to_identify.txt"

        with open(file_name, 'w') as fd:
            for results_name, segments_data in self._failed_to_identify.items():
                fd.write(results_name + "\n")
                for segment_name, generated_signatures in segments_data.items():
                    fd.write(f"\t\t{segment_name}:\t\t{generated_signatures}\n")
                fd.write("\n\n")


class StatsAggregatorManager:
    def __init__(self, simulation_folders: list[str], verbose: bool = True) -> None:
        """
        Initializes the StatsAggregator for further use.
        :param simulation_folders: list of simulation folders (containing results)
        :param verbose: if True, print progress information
        """

        self.simulation_folders: list[str] = simulation_folders
        self.verbose: bool = verbose

        self.result_files: dict[str, dict] = {}  # mapping of simulation folder name to results file contents (json dct)
        self._successfully_loaded: int = 0

        if self.verbose:
            print("[cyan][b]-> Loading workflow results...[/b][/cyan]", end="\r")

        self._load_results()

        if self.verbose:
            print("[green][b]✔  All available results were loaded              [/b][/green]")

    def _load_results(self):
        """
        Loads all results files from the simulation folders
        :return: None
        """
        for path in self.simulation_folders:
            results_path: str = os.path.join(path, "results.json")
            if not os.path.exists(results_path):
                continue
            with open(results_path, 'r') as fd:
                data_str: str = fd.read()
                if not data_str:
                    continue
                self.result_files[os.path.basename(path)] = json.loads(data_str)
            self._successfully_loaded += 1

    def run_aggregators(self, aggregator_kinds: list[str], verbose: bool = True) -> None:
        """
        Run aggregators given their names
        :param aggregator_kinds: types of aggregators to run (see StatsAggregatorKind)
        :param verbose: print messages about success/errors in parsing and running aggregators
        :return: None
        """
        aggregators_to_run: list[Type[StatsAggregator]] = parse_stats_types(aggregator_kinds, verbose)
        for agg in aggregators_to_run:
            agg(self).run()  # these are verbose by default


class StatsAggregatorKind(Enum):
    LOAD = LoadStatsAggregator
    SIGNATURES = SignatureStatsAggregator
    IDENTIFIERS = IdentifierStatsAggregator


def parse_stats_types(kinds: list[str], verbose: bool = True) -> list[Type[StatsAggregator]]:
    """
    Parses text-like aggregator kinds (case-insensitive) into aggregator classes
    :param kinds: to parse
    :param verbose: print messages about success/errors
    :return: list of TaskTypes
    """

    res: list[StatsAggregatorKind] = []
    for agg_txt in kinds:
        for agg_enum in StatsAggregatorKind:
            if agg_enum.name.lower() == agg_txt.lower():
                res.append(agg_enum)
                break
        else:
            if verbose:
                print(f"[red][b]✘  unable to parse aggregator '{agg_txt}', "
                      f"available options: {[x.name for x in StatsAggregatorKind]}[/b][/red]")
            return []


    if verbose:
        print(f"[green][b]✔  Following aggregators will be executed: "
              f"{[x.name for x in res]}[/b][/green]")
    return [x.value for x in res]