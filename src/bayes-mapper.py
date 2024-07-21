#!/usr/bin/env python3

from argparse import ArgumentParser, Namespace
import itertools
import os
from sys import argv, stderr, stdout
from typing import Any, Literal, Self
import typing
import Bio.Phylo  # type: ignore
from Bio.Phylo.BaseTree import Clade, Tree  # type: ignore

__version__ = "1.0.2"


TreeFormat = Literal["newick", "nexus", "phyloxml", "nexml"]
"""系統樹のフォーマットを表します。
"""


class Program:
    @classmethod
    def main(cls, args: list[str]) -> int:
        """main関数です。

        Args:
            args (list[str]): 引数

        Returns:
            int: Exit Code
        """
        parser: ArgumentParser = cls.create_parser()

        if len(args) == 0:
            parser.print_help(stdout)
            return 0

        arguments = Arguments(parser.parse_args(args))

        tree_format: TreeFormat = arguments.tree_format
        ml_tree: Tree = cls.read_tree(arguments.ml_tree, tree_format)
        bayes_tree: Tree = cls.read_tree(arguments.bayes_tree, tree_format)
        min_bp: int = arguments.min_bp
        min_pp: float = arguments.min_pp

        bayes_info = TreeInfo(bayes_tree)

        for _internal in ml_tree.get_nonterminals():
            ml_clade: Clade = _internal
            if ml_clade.confidence is None and ml_clade.name is None:
                continue
            bayes_clade: CladeInfo | None = bayes_info.find_same_clade(CladeInfo(ml_clade, ml_tree))
            bayes_confidence: str = "-" if bayes_clade is None or bayes_clade.confidence is None else str(bayes_clade.confidence)
            # フィルタリング
            # 条件に適合しない場合は評価を表示せず
            if (ml_clade.confidence is not None and ml_clade.confidence < min_bp) or (bayes_clade is not None and bayes_clade.confidence is not None and bayes_clade.confidence < min_pp):
                ml_clade.name = None
            else:
                ml_clade.name = f"{(ml_clade.name if ml_clade.confidence is None else str(ml_clade.confidence))}/{bayes_confidence}"
            ml_clade.confidence = None

        Bio.Phylo.write(ml_tree, arguments.out_path, tree_format, format_branch_length="%e")

        return 0

    @classmethod
    def create_parser(cls) -> ArgumentParser:
        """引数解析用のArgumentParserのインスタンスを生成します。

        Returns:
            ArgumentParser: ArgumentParserの新しいインスタンス
        """
        result = ArgumentParser()

        result.add_argument("-v", "--version", action="version", help="display the version", version=__version__)
        result.add_argument("-m", "--ml-tree", type=str, required=True, help="File of ML-Tree", metavar="FILE")
        result.add_argument("-b", "--bayes-tree", type=str, required=True, help="File of bayes tree", metavar="FILE")
        result.add_argument("-f", "--tree-format", type=str, required=False, default="newick", help="Format of tree file (default=newick)", metavar="STR")
        result.add_argument("-o", "--out", type=str, required=True, help="Output file name", metavar="FILE")
        result.add_argument("--min-bp", type=int, required=False, default=0, help="Minimum value of BP value", metavar="INT")
        result.add_argument("--min-pp", type=float, required=False, default=0, help="Minimum value of PP value", metavar="FLOAT")

        return result

    @classmethod
    def read_tree(cls, filename: str, format: str) -> Tree:
        result: Tree = Bio.Phylo.read(filename, format)
        cls.__assign_parent(result.root, None)
        return result

    @classmethod
    def __assign_parent(cls, target: Clade, parent: Clade | None) -> None:
        """再帰的に親ノードの情報を設定します。

        Args:
            target (Clade): 設定処理を行う対象クレード
            parent (Clade | None): 親クレード
        """
        target.parent: Clade | None = parent  # type: ignore
        for child in target.clades:
            cls.__assign_parent(child, target)


class TreeInfo:
    """ツリーの情報を表します。
    """

    def __init__(self, tree: Tree) -> None:
        """TreeInfoの新しいインスタンスを初期化します。

        Args:
            tree (Tree): BioPythonのTreeインスタンス
        """
        self.__clades: list[CladeInfo] = list[CladeInfo]()
        self.__init_clades(tree)

    def __init_clades(self, tree: Tree) -> None:
        """self.__cladesを初期化します。

        Args:
            tree (Tree): BioPythonのTreeのインスタンス
        """
        internals: list[Clade] = tree.get_nonterminals()
        for current in internals:
            self.__clades.append(CladeInfo(current, tree))

    def find_same_clade(self, target: "CladeInfo") -> "CladeInfo | None":
        """同値であるクレードを取得します。

        Returns:
            CladeInfo | None: 同値であるクレード。見つからない場合はNone
        """
        for current in self.__clades:
            if current.has_same_taxa(target):
                return current
        return None


class CladeInfo:
    """クレードの情報を表します。
    """

    @property
    def confidence(self) -> int | float | None:
        """信頼度を取得します。
        """
        return self.__confidence

    @property
    def name(self) -> str | None:
        """名前を取得します。
        """
        return self.__name

    def __init__(self, clade: Clade, tree: Tree) -> None:
        """CladeInfoの新しいインスタンスを初期化します。

        Args:
            clade (Clade): BioPythonのCladeのインスタンス
            tree (Tree): BioPythonのTreeのインスタンス
        """
        self.__confidence: int | float | None = clade.confidence
        self.__name: str | None = clade.name
        self.__taxa_list: list[list[str]] = list[list[str]]()
        self.__init_taxa_list(clade, tree)

    def __init_taxa_list(self, clade: Clade, tree: Tree) -> None:
        """self.__taxa_listを初期化します。

        Args:
            clade (Clade): BioPythonのCladeのインスタンス
            tree (Tree): BioPythonのTreeのインスタンス
        """
        children_names: list[str] = [self.__get_normalized_name(current) for current in clade.get_terminals()]
        children_names.sort()
        self.__taxa_list.append(children_names)
        children_set: set[str] = set(children_names)
        parent_set: set[str] = set([self.__get_normalized_name(current) for current in tree.get_terminals()]) - children_set
        parent_names: list[str] = list(parent_set)
        parent_names.sort()
        self.__taxa_list.append(parent_names)
        self.__taxa_list.sort()

    @staticmethod
    def __get_normalized_name(clade: Clade) -> str:
        clade_name: Any | None = clade.name
        if clade_name is None:
            return ""
        if isinstance(clade_name, str):
            return clade_name
        return str(clade_name)

    def has_same_taxa(self, target: Self) -> bool:
        """同値のTaxaを表すかどうかを検証します。

        Args:
            target (CladeInfo): 検証対象

        Remarks:
            このメソッドが処理の本体

            同じようなTaxaの隔て方をしているかどうかを見る。
            A, B, C, Dの系統樹を例とした場合，
            selfがA+BとC+Dに分ける時，targetが同様にA+BとC+Dに分けるかを検証している

        Returns:
            bool: targetとselfが同じ系統の分け方をしていたらTrue，それ以外でFalse
        """

        # 同じインスタンスなら常にTrue
        if self is target:
            return True

        return self.__taxa_list == target.__taxa_list

    def __str__(self) -> str:
        confidence: str = "-" if self.confidence is None else str(self.confidence)
        return f"{confidence}\t{str.join(' | ', [str.join(', ', current) for current in self.__taxa_list])}"


class ScriptAbortionError(Exception):
    """スクリプト内のエラーを表します。
    """
    @property
    def message(self) -> str:
        """エラーメッセージを取得します。
        """
        return self.__message

    @property
    def exit_code(self) -> int:
        """Exit Codeを取得します。
        """
        return self.__exit_code

    def __init__(self, message: str, exit_code: int | None = None):
        """ScriptAbortionErrorの新しいインスタンスを初期化します。

        Args:
            message (str): エラーメッセージ
            exit_code (int | None): Exit Code
        """
        super()
        self.__message: str = message
        self.__exit_code: int = 1 if exit_code is None else exit_code


class Arguments:
    """コマンド引数を表します。
    """

    @property
    def ml_tree(self) -> str:
        """MLツリーのパスを取得します。
        """
        result: str = self.__namespace.ml_tree

        if not os.path.isfile(result):
            raise ScriptAbortionError(f"ファイル'{result}'が存在しません")
        return result

    @property
    def bayes_tree(self) -> str:
        """ベイズツリーのパスを取得します。
        """
        result: str = self.__namespace.bayes_tree

        if not os.path.isfile(result):
            raise ScriptAbortionError(f"ファイル'{result}'が存在しません")
        return result

    @property
    def tree_format(self) -> TreeFormat:
        """ツリーフォーマットを取得します。
        """
        result: TreeFormat = self.__namespace.tree_format

        if result not in typing.get_args(TreeFormat):
            raise ScriptAbortionError(f"ツリーフォーマット'{result}'は無効です")
        return result

    @property
    def out_path(self) -> str:
        """出力先のパスを取得します。
        """
        return self.__namespace.out

    @property
    def min_bp(self) -> int:
        """フィルターするBP値の最小値を取得します。

        Returns:
            int: BP値の最小値
        """
        result: int = self.__namespace.min_bp
        if result < 0 or 100 < result:
            raise ScriptAbortionError(f"BP最小値'{result}'は0-100の整数値である必要があります")
        return result

    @property
    def min_pp(self) -> float:
        """フィルターするPP値の最小値を取得します。

        Returns:
            float: pP値の最小値
        """
        result: float = self.__namespace.min_pp
        if result < 0 or 1 < result:
            raise ScriptAbortionError(f"BP最小値'{result}'は0-1の値である必要があります")
        return result

    def __init__(self, namespace: Namespace) -> None:
        """Argumentsの新しいインスタンスを初期化します。

        Args:
            namespace (Namespace): 引数解析結果
        """
        self.__namespace: Namespace = namespace


if __name__ == "__main__":
    try:
        exit_code: int = Program.main(argv[1:])
        exit(exit_code)
    except ScriptAbortionError as e:
        print(e.message, file=stderr)
        exit(e.exit_code)
