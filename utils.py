from pathlib import Path
from datetime import datetime
from subprocess import run, STDOUT
from multiprocessing import Pool
import logging

# 创建日志记录器
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class Counter:

    def __init__(self):
        self.N = 0

    def step(self):
        self.N += 1
        return self.N


class Cmds:

    def __init__(self, cmds: list[str], shell_files: list[Path]):
        self.cmds = cmds
        self.shell_files = shell_files

    @staticmethod
    def done_or_dead(shell_file: Path, cmd: str):
        start_time = datetime.now()
        ok_file = shell_file.with_suffix(".ok")
        log_file = shell_file.with_suffix(".log")
        logger.info(f"{shell_file} STARTED!")
        try:
            if ok_file.exists():
                logger.info(f"{shell_file} SKIPPED!")
            else:
                with open(shell_file, "w") as fs:
                    fs.write(cmd)
                with open(log_file, "w") as fl:
                    result = run(["sh", str(shell_file)], stdout=fl, stderr=STDOUT)
                end_time = datetime.now()
                if result.returncode == 0:
                    ok_file.touch()
                    logger.info(f"{shell_file} DONE! ELAPSED: {end_time - start_time}s.")
                else:
                    logger.error(f"{shell_file} FAILED! ELAPSED: {end_time - start_time}s.")
                    logger.error(result.stderr.decode())
                    return False
            return True
        except Exception as ex:
            logger.error(f"Exception occurred: {ex}")
            return False

    def multi_run(self, cpu) -> None:

        with Pool(cpu) as pool:
            results = [
                pool.apply_async(self.done_or_dead, (shell_file, cmd))
                for shell_file, cmd in zip(self.shell_files, self.cmds)
            ]
            pool.close()
            pool.join()

        failed_files = [shell_file for shell_file, result in zip(self.shell_files, results) if not result.get()]
        if failed_files:
            raise Exception(f"The following files failed: {failed_files}")


class Results:

    def __init__(self, files: list[Path]):
        self.files = files

    def check_empty(self):
        for file in self.files:
            if file.stat().st_size == 0:
                raise Exception(f"Error: FILE:{file} IS EMPTY!")
        return self

    def check_exists(self):
        for file in self.files:
            if not file.exists():
                raise Exception(f"Error: FILE:{file} DOES NOT EXIST!")
        return self
