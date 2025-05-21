import subprocess

def run_conda_command(args):
    try:
        subprocess.run(["conda"] + args, check=True)
        print(f"✅ Successfully ran: {' '.join(args)}")
    except subprocess.CalledProcessError:
        print(f"❌ Error running: {' '.join(args)}")

def main():
    env_name = "comolaenv"

    # 1. Create conda environment with Python 3.11
    run_conda_command(["create", "-n", env_name, "python=3.11", "-y"])

    # 2. Install dependencies in the new environment
    dependencies = [
        ("r-base", "conda-forge"),
        ("r-essentials", "conda-forge"),
        ("matplotlib", "conda-forge"),
        ("numpy=1.26.4", "conda-forge"),
    ]

    for package, channel in dependencies:
        run_conda_command(["install", "-n", env_name, "-c", channel, package, "-y"])

    print(f"\n🎉 Environment '{env_name}' is ready! To activate it, run:\n\n    conda activate {env_name}\n")
    print("🔍 Then you can check paths using:\n\n    where python\n    where R\n")

if __name__ == "__main__":
    main()
