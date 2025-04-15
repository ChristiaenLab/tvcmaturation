{

  nixConfig = {
    bash-prompt = "\[TVCmaturation$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

	BSgenomeKY = {
	  url = "github:kewiechecki/BSgenome.Crobusta.HT.KY";
	  flake = false;
	};
	CrobustaTFs = {
	  url = "github:kewiechecki/CrobustaTFs";
	  flake = false;
	};
	dirfns = {
	  url = "github:kewiechecki/dirfns";
	  flake = false;
	};
	moreComplexHeatmap = {
	  url = "github:kewiechecki/moreComplexHeatmap";
	  flake = false;
	};
	peakToGene = {
	  url = "github:kewiechecki/peakToGene";
	  flake = false;
	};
	tfenrichr = {
	  url = "github:kewiechecki/tfenrichr";
	  flake = false;
	};
  };

  outputs = { self, nixpkgs, utils, REPLVim, igraph_jll, leiden_jll, Leiden,
              Autoencoders, TrainingIO, DictMap, DeePWAK }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          #config.cudaSupport = system == "x86_64-linux";
        };

		rScript = ''
library(devtools)

install_github("__BSGENOMEKY__")
install_github("__CROBUSTATFS__")
install_github("__DIRFNS__")
install_github("__MORECOMPLEXHEATMAPS__")
install_github("__PEAKTOGENE__")
install_github("__TFENRICHR__")
'';
      in {
        defaultPackage = pkgs.stdenv.mkDerivation {
          name = "CrobustaScreen";
          src = ./.;
          buildInputs = [];
        };

        devShell = with pkgs; mkShell {
          name = "crobusta-screen-shell";
          buildInputs = [
            R
            pkgs.rPackages.devtools
			pkgs.rPackages.optparse
            pkgs.rPackages.purrr

            pkgs.rPackages.GenomicFeature
            pkgs.rPackages.motifmatchr
            pkgs.rPackages.TFBSTools
            
			# to fetch interactions
            pkgs.rPackages.biomaRt
            pkgs.rPackages.STRINGdb
			# clustering
            #pkgs.rPackages.class
            #pkgs.rPackages.cluster
            #pkgs.rPackages.fgsea
            #pkgs.rPackages.igraph
            #pkgs.rPackages.leiden
			# visualization
			pkgs.rPackages.circlize
			pkgs.rPackages.ComplexHeatmap
			pkgs.rPackages.ggplot2
			pkgs.rPackages.ggpubr
			#pkgs.rPackages.umap
            # Python environment with selected packages.
            #(python3.withPackages (ps: with ps; [ umap-learn leidenalg igraph ]))
            #python3Packages.virtualenv
            #julia
            git      # for git prompt support
            #pkgs.stdenv.cc
            #pkgs.gfortran
          ];
          shellHook = "
source ${git}/share/bash-completion/completions/git-prompt.sh

cat > deps.R <<'EOF'
${rScript}
EOF

sed -i \"s|__BSGENOMEKY__|${toString BSgenomeKY}|g\" deps.R
sed -i \"s|__CROBUSTATFS__|${toString CrobustaTFs}|g\" deps.R
sed -i \"s|__DIRFNS__|${toString dirfns}|g\" deps.R
sed -i \"s|__MORECOMPLEXHEATMAPS__|${toString moreComplexHeatmaps}|g\" deps.R
sed -i \"s|__PEAKTOGENE__|${toString peakToGene}|g\" deps.R
sed -i \"s|__TFENRICHR__|${toString tfenrichr}|g\" deps.R

# Replace the dollar placeholder with a literal dollar sign.
sed -i \"s|__DOLLAR_PLACEHOLDER__|\\\$|g\" deps.R

Rscript deps.R
";
		};
	});
}
