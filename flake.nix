{
  description = "Cpp Environment";

  nixConfig = {
    bash-prompt-prefix = "(CppEnv) ";
    bash-prompt = ''\[\033]0;\u@\h:\w\007\]\[\033[01;32m\]\u@\h\[\033[01;34m\] \w \$\[\033[00m\]'';
    bash-prompt-suffix = " ";
  };
  
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    devenv.url = "github:cachix/devenv/latest";
  };

  
  outputs = { self, nixpkgs, flake-utils, devenv, ... } @ inputs:
    flake-utils.lib.eachSystem [
      "x86_64-linux"
    ]
      (system:
        let

          nixpkgs-patched = (import nixpkgs { inherit system; }).applyPatches {
            name = "nixpkgs-patched";
            src = nixpkgs;
            patches = [ ./patches/nix_fix_linkerscript.diff ];
          };

          pkgs = import nixpkgs-patched { inherit system; };

      in rec {
        devShell = devenv.lib.mkShell {
          inherit inputs pkgs;
          
          modules = [{
    
            scripts.bazel.exec = (''

              BAZEL_CXXOPTS="-isystem${pkgs.ginac}/include:-isystem${pkgs.arb}/include:-isystem${pkgs.flint}/include:-isystem${pkgs.cln}/include:-isystem${pkgs.boost.dev}/include:-isystem${pkgs.gmp.dev}/include:-isystem${pkgs.mpfr.dev}/include"
              BAZEL_LINKOPTS="-L${pkgs.ginac}/lib:-L${pkgs.cln}/lib:-L${pkgs.flint}/lib:-L${pkgs.arb}/lib"

              if [[
                "$1" == "build" ||
                "$1" == "run"
              ]]; then 
                bazelisk $1 \
                  --action_env=BAZEL_CXXOPTS=$BAZEL_CXXOPTS \
                  --action_env=BAZEL_LINKOPTS=$BAZEL_LINKOPTS \
                  --experimental_enable_bzlmod \
                  ''${@:2}
              else
                bazelisk $@
              fi
            '');
            
            
            packages = [
              pkgs.arb 
              pkgs.flint
              pkgs.ginac
              pkgs.gmp
              pkgs.cln
              pkgs.mpfr
              pkgs.gcc12
              pkgs.gcc12Stdenv
              pkgs.bazelisk
              pkgs.boost
            ]; 

            enterShell = ''
              # Ensure that the bazel command points to our custom wrapper.
              [[ $(type -t bazel) == "alias" ]] && unalias bazel


              # Prevent rules_cc from using anything other than clang.
              # export CC=clang

              # Probably a bug in nix. Setting LD=ld.lld here doesn't work.
              # export LD=${pkgs.llvmPackages_15.lld}/bin/ld.lld

              # Prettier color output for the ls command.
              alias ls='ls --color=auto'

              echo ${pkgs.gcc12}
            '';
          }];
        };
    });
}
