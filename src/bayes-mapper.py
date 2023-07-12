from sys import argv


class Program:
    @classmethod
    def main(cls, args: list[str]) -> int:
        """main関数です。

        Args:
            args (list[str]): 引数

        Returns:
            int: Exit Code
        """
        return 0


if __name__ == "__main__":
    exit_code: int = Program.main(argv[1:])
    exit(exit_code)
