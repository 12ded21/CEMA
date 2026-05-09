import os
import subprocess
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_EXECUTABLE = PROJECT_ROOT / "build" / "main"
DEFAULT_DATA_DIR = PROJECT_ROOT / "data" / "DIMACS-10"

def resolve_project_path(path):
    path = Path(path)
    return path if path.is_absolute() else PROJECT_ROOT / path

def run_task(executable, time_limit, d, m_val, e_val, file_path, output_file, file_name):
    print(f"\nProcessing file: {file_name}")
    print(f"  Running with d={d}, m={m_val}, edge_method={e_val}, appending to {output_file}")
    
    # 构建命令时，把 -m 和 -e 参数都加进去
    cmd = [
        executable,
        "-l", str(time_limit),
        "-d", str(d),
        "-m", str(m_val),
        "-e", str(e_val),
        str(file_path),
        str(output_file)
    ]
    
    # 执行命令
    try:
        subprocess.run(cmd, check=True)
        print(f"  Successfully completed d={d}, m={m_val}, e={e_val} for {file_name}")
    except subprocess.CalledProcessError as e:
        print(f"  Error running command for d={d}, m={m_val}, e={e_val}: {e}") 
    except Exception as e:
        print(f"  Unexpected error for d={d}, m={m_val}, e={e_val}: {e}")

def run_experiments(data_dirs, executable, d_values, m_values, e_values, time_limit, threads):
    # 提前为每个 d, m_val 和 e_val 的组合创建/清空结果文件
    for d in d_values:
        for m_val in m_values:
            for e_val in e_values:
                output_file = f"d-sweep-m={m_val}-d={d}-edge_method={e_val}.txt"
                with open(output_file, 'w') as f:
                    pass  # 清空文件内容
                print(f"Initialized output file: {output_file}")

    # 用来收集所有需要并行的任务
    tasks = []

    for data_dir in data_dirs:
        data_path = Path(data_dir)
        if not data_path.exists():
            print(f"Warning: data directory does not exist: {data_path}")
            continue
        all_files = list(data_path.rglob('*'))
        print(f"Found {len(all_files)} files in {data_dir}")
        
        for file_path in all_files:
            if not file_path.is_file():
                continue
                
            # 👉 严密过滤：跳过所有 .txt, .log, .md 结尾的文件，或以 d-sweep- 开头的文件
            # 彻底杜绝 C++ 因为读取非图文件而陷入死循环！
            if file_path.suffix in ['.txt', '.log', '.md'] or file_path.name.startswith("d-sweep-"):
                continue
                
            file_name = file_path.name
            
            # 三重循环：对每个 d 值, m 值和 e 值收集任务
            for d in d_values:
                for m_val in m_values:
                    for e_val in e_values:
                        output_file = f"d-sweep-m={m_val}-d={d}-edge_method={e_val}.txt"
                        tasks.append((executable, time_limit, d, m_val, e_val, file_path, output_file, file_name))

    print(f"\nStarting {len(tasks)} tasks with {threads} threads...")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for task_args in tasks:
            executor.submit(run_task, *task_args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run batch experiments by sweeping budget density d')
    parser.add_argument('--executable', default=str(DEFAULT_EXECUTABLE), help='Path to the executable')
    parser.add_argument('--data-dir', nargs='+', default=[DEFAULT_DATA_DIR], help='Data directory or directories, relative paths are resolved from the project root')
    parser.add_argument('--time-limit', type=int, default=605, help='Time limit in seconds')
    # 默认线程数设置为了 96
    parser.add_argument('--threads', type=int, default=16, help='Number of concurrent threads')
    # 接受多个 m 方法的列表，默认执行 1, 2, 3, 4, 5
    parser.add_argument('--m_methods', nargs='+', type=int, default=[2], help='List of LB selection methods (-m)')
    # 接受多个 e 方法的列表，默认执行 5
    parser.add_argument('--edge_methods', nargs='+', type=int, default=[5], help='List of edge selection methods (-e)')
    
    args = parser.parse_args()
    
    data_dirs = [resolve_project_path(path) for path in args.data_dir]
    
    # 你的 d_values 列表
    d_values = [0.0001, 0.0005, 0.001, 0.0015, 0.002]
    
    # 确保传入的 m_values 都是合法的 (1-5 之间)
    m_values = [m for m in args.m_methods if 1 <= m <= 5]
    if not m_values:
        print("Error: No valid m methods provided. Please use numbers between 1 and 5.")
        exit(1)

    # 确保传入的 e_values 都是合法的 (1-5 之间)
    e_values = [e for e in args.edge_methods if 1 <= e <= 5]
    if not e_values:
        print("Error: No valid edge methods provided. Please use numbers between 1 and 5.")
        exit(1)
    
    print("Starting experiments...")
    print(f"Executable: {args.executable}")
    print(f"Time limit: {args.time_limit} seconds")
    print(f"Threads: {args.threads}")
    print(f"d values: {d_values}")
    print(f"m values (LB Methods): {m_values}")
    print(f"e values (Edge Methods): {e_values}")
    
    # 调用时把 d_values, m_values 和 e_values 传进去
    run_experiments(data_dirs, args.executable, d_values, m_values, e_values, args.time_limit, args.threads)
    print("\nAll experiments completed!")
