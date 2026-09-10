\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Liebert_Author_revised}[2026/09/08 v1.01
^^J LaTeX class file for Liebert Journals applicable for preparing manuscript at author's end]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[12pt,twoside]{article}

\usepackage{mathtools,amsthm,amsfonts,amssymb}%
\usepackage[round,sort]{natbib}%
\usepackage[doublespacing]{setspace}%
\usepackage{threeparttable,booktabs}%
\usepackage{graphicx}%
\usepackage[figuresright]{rotating}%
\usepackage{enumitem}%

\newcommand{\keywords}[1]{\bgroup\noindent\textbf{Keywords:}
#1\egroup}%

\renewenvironment{abstract}{%
\begin{quote}\textbf{Abstract}\par\noindent }%
{\end{quote}}%

\usepackage{hyperref}%

\usepackage[ruled,norelsize]{algorithm2e}%

\usepackage{mathrsfs}%

\newtheorem{theorem}{Theorem}%
\newtheorem{lemma}{Lemma}%
\newtheorem{corollary}{Corollary}%
\newtheorem{remark}{Remark}%

\newtheorem*{theorem*}{Theorem}%
\newtheorem*{lemma*}{Lemma}%
\newtheorem*{corollary*}{Corollary}%
\newtheorem*{remark*}{Remark}%

\AtBeginDocument{\pagestyle{headings}}%

\newcommand{\notetocomp}[1]{\bgroup\parindent\z@\raggedright%
\marginpar{\footnotesize#1\par}
\egroup}%

\endinput
