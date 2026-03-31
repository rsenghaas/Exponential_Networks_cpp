{
  description = "Cpp Environment";

  nixConfig = {
    bash-prompt-prefix = "(CppEnv) ";
    bash-prompt =
      ''\[\033]0;\u@\h:\w\007\]\[\033[01;32m\]\u@\h\[\033[01;34m\] \w \$\[\033[00m\]'';
    bash-prompt-suffix = " ";
  };

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    devenv.url = "github:cachix/devenv";
    pre-commit-hooks-nix.url = "github:cachix/pre-commit-hooks.nix";
  };

  outputs = inputs@{ self, nixpkgs, flake-utils, devenv, pre-commit-hooks-nix, ... }:
    flake-utils.lib.eachSystem [ "x86_64-linux" ] (system:
      let
        system = "x86_64-linux";
        pkgs = import nixpkgs { inherit system; };

        hooks = import ./pre-commit-hooks.nix {
          inherit pkgs;
        };

        # ✔️ Proper nixpkgs python env (no parsing hacks)
        pythonEnv = pkgs.python3.withPackages (ps: with ps; [
          matplotlib
          numpy
          pandas
          sympy
          networkx
          scipy
          pypdf
          # add more packages here explicitly
        ]);

      in {
        checks.pre-commit-check =
          pre-commit-hooks-nix.lib.${system}.run {
            src = ./.;
            hooks = hooks;
          };

        devShells.default = devenv.lib.mkShell {
          inherit inputs pkgs;

          modules = [
            ({ pkgs, config, ... }: { 
            languages.python.enable = true;
            languages.python.package = pythonEnv;


            scripts.bazelnr.exec = ''
              echo "In bazel non recursive"

              BAZEL_CXXOPTS="-isystem${pkgs.ginac}/include:-isystem${pkgs.flint}/include:-isystem${pkgs.cln}/include:-isystem${pkgs.boost.dev}/include:-isystem${pkgs.gmp.dev}/include:-isystem${pkgs.mpfr.dev}/include"
              BAZEL_LINKOPTS="-L${pkgs.ginac}/lib:-L${pkgs.cln}/lib:-L${pkgs.flint}/lib"

              if [[ "$1" == "build" || "$1" == "run" ]]; then
                bazel $1 \
                  --action_env=BAZEL_CXXOPTS=$BAZEL_CXXOPTS \
                  --action_env=BAZEL_LINKOPTS=$BAZEL_LINKOPTS \
                  --experimental_enable_bzlmod \
                  ''${@:2}
              else
                bazel $@
              fi
            '';

            packages = with pkgs; [
              # arb
              fmt
              ripgrep
              flint
              ginac
              gmp
              cln
              spdlog
              mpfr
              stdenv.cc
              bazel
              boost
              zsh
            ];

            enterShell = ''
              export SHELL=${pkgs.zsh}/bin/zsh
              exec ${pkgs.zsh}/bin/zsh
              alias ls='ls --color=auto'
            '';
            })
          ];
        };
      });
}

