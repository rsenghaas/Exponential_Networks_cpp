{
  description = "Cpp Environment";

  nixConfig = {
    bash-prompt-prefix = "(CppEnv) ";
    bash-prompt = ''\[\033]0;\u@\h:\w\007\]\[\033[01;32m\]\u@\h\[\033[01;34m\] \w \$\[\033[00m\]'';
    bash-prompt-suffix = " ";
  };
  
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    devenv.url = "github:cachix/devenv/latest";
    pre-commit-hooks-nix.url = "github:cachix/pre-commit-hooks.nix";
    mach-nix.url = "github:davhau/mach-nix";
  };

  
  outputs = { self
            , nixpkgs
            , flake-utils
            , devenv
            , mach-nix
            , pre-commit-hooks-nix
            , ... } @ inputs:
    flake-utils.lib.eachSystem [
      "x86_64-linux"
    ]
      (system:
        let
          mach = mach-nix.lib.${system};
          pythonEnv = mach.mkPython {
            python = "python310";
            requirements = builtins.readFile ./requirements.txt;
          };
          nixpkgs-patched = (import nixpkgs { inherit system; }).applyPatches {
            name = "nixpkgs-patched";
            src = nixpkgs;
            patches = [ ./patches/nix_fix_linkerscript.diff ];
          };
          
          pkgs = import nixpkgs { inherit system; };
      in rec {
        hooks = import ./pre-commit-hooks.nix {
          inherit pkgs;
        };

        checks = {
          pre-commit-check = pre-commit-hooks-nix.lib.${system}.run {
            src = ./.;
            hooks = hooks;
          };
        };

        devShell = devenv.lib.mkShell {
          inherit inputs pkgs;
          modules = [{
          scripts.bazelnr.exec = (''
              echo "In bazel non recursive"
              BAZEL_CXXOPTS="-isystem${pkgs.ginac}/include:-isystem${pkgs.arb}/include:-isystem${pkgs.flint}/include:-isystem${pkgs.cln}/include:-isystem${pkgs.boost.dev}/include:-isystem${pkgs.gmp.dev}/include:-isystem${pkgs.mpfr.dev}/include"
              BAZEL_LINKOPTS="-L${pkgs.ginac}/lib:-L${pkgs.cln}/lib:-L${pkgs.flint}/lib:-L${pkgs.arb}/lib"
              if [[
                "$1" == "build" ||
                "$1" == "run"
              ]]; then 
                bazel $1 \
                  --action_env=BAZEL_CXXOPTS=$BAZEL_CXXOPTS \
                  --action_env=BAZEL_LINKOPTS=$BAZEL_LINKOPTS \
                  --experimental_enable_bzlmod \
                  ''${@:2}
              else
                bazel $@
              fi
            '');
            
            packages = [
              pkgs.arb 
              pkgs.fmt
              pkgs.ripgrep
              pkgs.flint
              pkgs.ginac
              pkgs.gmp
              pkgs.cln
              pkgs.spdlog
              pkgs.mpfr
              pkgs.gcc13
              pkgs.gcc13Stdenv
              # pkgs.bazelisk
              pkgs.bazel
              pkgs.boost
              pythonEnv
            ]; 

            enterShell = ''
              alias ls='ls --color=auto'
              export PYTHONPATH="${pythonEnv}/bin/python"
            '';
          }];
        };
    });
}
