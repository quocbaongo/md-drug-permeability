# md-drug-permeability
\documentclass[11pt, a4paper]{article}

% Essential packages for scientific writing
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath, amssymb}  % For advanced math and symbols
\usepackage{geometry}           % For page margins
\usepackage{hyperref}           % For hyperlinks
\usepackage{enumitem}           % For better list formatting
\usepackage{xcolor}             % For coloring links (optional)

% Page setup
\geometry{margin=1in}

% Hyperlink setup
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,      
    urlcolor=blue,
    pdftitle={MD Permeability Workflow},
    pdfpagemode=FullScreen,
}

\title{\textbf{Permeability Coefficient Estimation Workflow: \\ Drug-Membrane Interactions}}
\author{Computational Biology Research Group}
\date{\today}

\begin{document}

\maketitle

\section{Overview}
This document outlines the computational workflow for theoretically estimating the permeability coefficient of small molecule drugs across lipid bilayers using Molecular Dynamics (MD) simulations. The pipeline focuses on the researched drug molecule \textbf{YT0} (\href{https://www.rcsb.org/ligand/YT0}{RCSB Ligand ID: YT0}) and implements the Inhomogeneous Solubility-Diffusion Model (ISDM).

The workflow is divided into three primary modules:
\begin{enumerate}
    \item \textbf{Parameterization:} Derivation of high-fidelity simulation parameters.
    \item \textbf{Validation:} Partition coefficient ($\log P_{ow}$) calculation.
    \item \textbf{Permeability:} Computation of Free Energy Profiles (FEP) and Local Diffusion Coefficients to determine membrane resistance.
\end{enumerate}

\section{Software \& Dependencies}
The workflow relies on the following software stack:
\begin{itemize}
    \item \textbf{MD Engine:} GROMACS (v2022.2)
    \item \textbf{Parameterization:} AmberTools2024, Gaussian16 (for QM calculations)
    \item \textbf{Analysis \& Scripting:} Python 3.x
    \begin{itemize}
        \item \textit{Libraries:} \texttt{numpy}, \texttt{rdkit}, \texttt{cclib}, \texttt{matplotlib}
    \end{itemize}
\end{itemize}

\section{Methodology}

\subsection{1. Drug Molecule Parameterization}
\textbf{Script:} \texttt{Parameterization-workflow.sh}

This module provides a detailed protocol for generating topology and coordinate files for the drug molecule. The parameterization strategy ensures high accuracy for non-standard residues:
\begin{itemize}
    \item \textbf{Partial Charges:} Calculated via Quantum Mechanical (QM) methods using Gaussian16. Charges are typically fitted using the RESP (Restrained Electrostatic Potential) or equivalent scheme.
    \item \textbf{Bonded Parameters:} Derived from the General Amber Force Field (GAFF).
\end{itemize}

\subsection{2. Validation: Water-Octanol Partitioning ($\log P_{ow}$)}
To validate the force field parameters derived in step 1, we compute the octanol-water partition coefficient ($\log P_{ow}$). This is achieved by calculating the solvation free energies of the drug in water ($\Delta G_{\text{solv}}^{\text{water}}$) and in octanol ($\Delta G_{\text{solv}}^{\text{oct}}$).

The partition coefficient is calculated as:
\begin{equation}
    \log P_{ow} = \frac{\Delta G_{\text{solv}}^{\text{water}} - \Delta G_{\text{solv}}^{\text{oct}}}{RT \ln(10)}
\end{equation}
Where:
\begin{itemize}
    \item $R$ is the universal gas constant ($8.314 \times 10^{-3} \, \text{kJ} \cdot \text{mol}^{-1} \cdot \text{K}^{-1}$).
    \item $T$ is the absolute temperature.
    \item $\ln(10) \approx 2.30258$ conversion factor.
\end{itemize}

\subsection{3. Permeability Coefficient Computation}
\textbf{Reference:} \href{https://www.sciencedirect.com/science/article/pii/S0006349522007378}{Biophysical Journal, 2022 (DOI: 10.1016/j.bpj.2022.08.025)}

This module implements the workflow for calculating the membrane permeability coefficient ($P$). Using the Inhomogeneous Solubility-Diffusion Model, the total resistance to permeation ($R$) is obtained by integrating the local resistance across the bilayer normal ($z$).

The relationship is defined as:
\begin{equation}
    \frac{1}{P} = R = R_{\text{hydrodynamic}} + \int_{z_{1}}^{z_{2}} \frac{e^{\beta \Delta G(z)}}{D(z)} \, dz
\end{equation}
Where:
\begin{itemize}
    \item $R_{\text{hydrodynamic}}$: Resistance contribution from the hydrodynamic water layer (approximated as $30 \pm 6$ in specific setups).
    \item $\Delta G(z)$: Potential of Mean Force (PMF) relative to bulk water.
    \item $D(z)$: Local diffusion coefficient profile along the membrane normal.
    \item $\beta = (k_B T)^{-1}$.
\end{itemize}

\end{document}
