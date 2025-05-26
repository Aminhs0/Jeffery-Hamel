\documentclass[12pt]{article}
\usepackage{multirow}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{hyperref}
\usepackage{subcaption}

\begin{document}
\fontsize{12pt}{14pt}\selectfont

% Assuming etc/head.tex includes necessary configurations
\input{etc/head}

\begin{abstract}
\noindent
This paper presents a Python-based numerical analysis of Jeffery-Hamel flow, a classical problem in fluid mechanics. The study considers a steady, radial flow of an incompressible Newtonian fluid. In the first part, the Navier-Stokes equations and boundary conditions are simplified to analyze the dimensionless velocity profile, pressure distribution along the centerline, and stress components on the channel walls for various positive and negative Reynolds numbers. In the second part, the energy conservation equation is simplified using a temperature function to derive the governing equation for the temperature distribution in the flow. Results are visualized using Python's Matplotlib library, providing insights into the flow and thermal behavior under different conditions.
\end{abstract}

\section{Introduction}
Jeffery-Hamel flow describes the motion of a viscous fluid between two converging or diverging planar walls at a given angle. This problem is fundamental in fluid mechanics due to its rich mathematical structure and applications in lubrication theory, aerodynamics, and microfluidics.

While early studies used simplifying assumptions or perturbation methods for approximate solutions, advances in computational techniques have enabled more accurate analysis. Most prior work focuses on flow characteristics, with fewer studies on thermal effects. This paper simplifies the governing equations for Jeffery-Hamel flow and implements them in Python, analyzing velocity, pressure, stress, and temperature distributions.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.4\textwidth]{figs/image.png}
\caption{Configuration of the Jeffery-Hamel flow channel.}
\label{fig:jeffery_hamel_config}
\vspace{10pt}
\end{figure}

\section{Governing Equations}
The Jeffery-Hamel flow is governed by mass, momentum, and energy conservation. Assuming steady, incompressible, axisymmetric flow, and using a similarity variable, the equations reduce to a system of ordinary differential equations (ODEs) in radial coordinates.

The continuity equation is given by:
\begin{equation} \label{eq:continuity}
\frac{1}{r} \frac{\partial}{\partial r} (r v_r) = 0
\end{equation}
which implies:
\begin{equation} \label{eq:velocity_function}
r v_r = f(\theta)
\end{equation}
At the centerline (\(\theta=0\)):
\begin{equation} \label{eq:centerline_velocity}
f(0) = C_0 = v_0 r
\end{equation}
with \(v_0\) and \(C_0\) positive for outflow, negative for inflow. The dimensionless velocity is defined as:
\begin{equation} \label{eq:dimensionless_velocity}
v^* = \frac{v_r}{v_0} = F(\theta)
\end{equation}

The momentum equation in the \(\theta\)-direction:
\begin{equation} \label{eq:momentum_theta}
\frac{\partial p}{\partial \theta} = \frac{2\mu}{r} \frac{\partial v_r}{\partial \theta}
\end{equation}
and in the \(r\)-direction:
\begin{equation} \label{eq:momentum_r}
\rho v_r \frac{\partial v_r}{\partial r} = -\frac{\partial p}{\partial r} + \mu \left[ \frac{1}{r} \frac{\partial v_r}{\partial r} + \frac{\partial^2 v_r}{\partial r^2} - \frac{v_r}{r^2} \right]
\end{equation}
Using (\ref{eq:dimensionless_velocity}) and (\ref{eq:momentum_theta}) leads to:
\begin{equation} \label{eq:velocity_derivative}
\frac{\left[ -\frac{1}{r} \left( \frac{\partial v_r}{\partial \theta} \right) \right]}{\frac{v_0}{r \alpha}} = -F'(\eta)
\end{equation}
where \(\eta = \theta/\alpha\). Substituting and eliminating pressure yields a third-order nonlinear ODE:
\begin{equation} \label{eq:jeffery_hamel_ode}
F'''(\eta) + (2\alpha \mathrm{Re}) F(\eta) F'(\eta) + (4 \alpha^2) F'(\eta) = 0
\end{equation}
where:
\begin{itemize}
    \item \(F\): Dimensionless velocity
    \item \(\alpha\): Half-angle between the plates
    \item \(\eta\): Dimensionless radial coordinate
    \item \(\mathrm{Re}\): Reynolds number
\end{itemize}

The pressure, shear stress, and normal stress coefficients are given by:
\begin{equation} \label{eq:pressure_integral}
p_\infty - p = \int_{\infty,0}^{r,0} \frac{\partial p}{\partial r} dr + \int_{r,0}^{r,\eta} \frac{\partial p}{\partial \eta} d\eta
\end{equation}
\begin{equation} \label{eq:pressure_coefficient}
c_p = \frac{p_\infty - p}{\frac{1}{2} \rho v_0^2} = 1 + \frac{4\alpha^2}{\alpha \mathrm{Re}} \left[ 1 - F(\eta; \alpha, \mathrm{Re}) \right] + \frac{1}{\alpha \mathrm{Re}} F'''(0; \alpha, \mathrm{Re})
\end{equation}
\begin{equation} \label{eq:shear_stress}
c_f = \frac{\tau_{\theta r}}{\frac{1}{2} \rho v_0^2} = \frac{2}{\mathrm{Re}} F'(\eta)
\end{equation}
\begin{equation} \label{eq:normal_stress}
c_n = \frac{\tau_{rr}}{\frac{1}{2} \rho v_0^2} = -\frac{\tau_{\theta\theta}}{\frac{1}{2} \rho v_0^2} = \frac{4\alpha}{\mathrm{Re}} F(\eta)
\end{equation}

The energy equation is:
\begin{equation} \label{eq:energy}
\rho C_p \left( v_r \frac{\partial T}{\partial r} \right) = k \left[ \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial T}{\partial r} \right) + \frac{1}{r^2} \frac{\partial^2 T}{\partial \theta^2} \right] + \Phi
\end{equation}
where \(C_p\) is the specific heat, \(k\) the thermal conductivity, and \(\Phi\) the viscous dissipation.

Defining a temperature function:
\begin{equation} \label{eq:temperature_function}
T - T_s = \frac{v_0^2 G(\theta)}{C_p r^2}
\end{equation}
and substituting into the energy equation, we obtain:
\begin{align}
G'' + (4 + 2 F \mathrm{Pr}) G + \mathrm{Pr} (4 F^2 - F'^2) = 0 \label{eq:temperature_ode_with_dissipation} \\
G'' + (4 + 2 F \mathrm{Pr}) G = 0 \label{eq:temperature_ode_no_dissipation}
\end{align}
where \(\mathrm{Pr}\) is the Prandtl number.

\subsection{Boundary Conditions}
Boundary conditions are specified at \(\eta=0\) and \(\eta=1\). At \(\eta=0\), the velocity matches the free stream value and the temperature equals the wall value, enforcing the no-slip and constant wall temperature conditions:
\begin{equation}
F(0) = 0, \quad F'(0) = 1, \quad G(0) = 0
\end{equation}
At \(\eta=1\):
\begin{equation}
F(1) = 0, \quad G(\alpha) = 0
\end{equation}

Parameters and boundary conditions are summarized in Table~\ref{tab:parameters}.
\begin{table}[!htbp]
\centering
\begin{tabular}{|c|c|}
\hline
Parameter & Value \\
\hline
\(F'(0)\) & 1 \\
\(F(0)\) & 0 \\
\(F(1)\) & 0 \\
\(2\alpha\) & 20$^\circ$ \\
\(G'(0)\) & 0 \\
\(G(\alpha)\) & 0 \\
\(\rho\) & 1.225 kg/m$^3$ \\
\(v_0\) & 1 m/s \\
\(C_p\) & 1005 J/kg$\cdot$K \\
\(k\) & 0.0262 W/m$\cdot$K \\
\(p_\infty\) & 101325 Pa \\
\(T_s\) & 300 K \\
\hline
\end{tabular}
\caption{Parameters and boundary conditions for Jeffery-Hamel flow.}
\label{tab:parameters}
\vspace{10pt}
\end{table}

\section{Method and Solution Approach}
The Jeffery-Hamel flow problem is numerically solved as a boundary value problem (BVP) using Python's \texttt{scipy.integrate.solve\_bvp}, which uses a fourth-order Runge-Kutta method with automatic mesh selection. The ODEs in Equations~\ref{eq:jeffery_hamel_ode}, \ref{eq:temperature_ode_with_dissipation}, and \ref{eq:temperature_ode_no_dissipation} are solved.

\subsection{Numerical Method}
The \texttt{solve\_bvp} function is designed for first-order BVPs; thus, the third-order ODE (\ref{eq:jeffery_hamel_ode}) is rewritten as a system of first-order ODEs. The Runge-Kutta scheme is:
\begin{align*}
k_1 &= h f(t_n, y_n), \\
k_2 &= h f(t_n + \frac{h}{2}, y_n + \frac{k_1}{2}), \\
k_3 &= h f(t_n + \frac{h}{2}, y_n + \frac{k_2}{2}), \\
k_4 &= h f(t_n + h, y_n + k_3), \\
y_{n+1} &= y_n + \frac{1}{6} (k_1 + 2k_2 + 2k_3 + k_4).
\end{align*}

\subsection{Parameter Definition}
Parameters include the channel half-angle \(\alpha = 10^\circ\) (in radians), Reynolds numbers (\(\mathrm{Re} = [100, 50, 5, -5, -50, -100]\)), density (\(\rho = 1.225\, \mathrm{kg/m^3}\)), specific heat (\(C_p = 1005\, \mathrm{J/kg\,K}\)), and other constants given in Table~\ref{tab:parameters}.

\subsection{ODE and Boundary Conditions}
The Jeffery-Hamel ODE is defined via a function, taking \(\eta\), the dependent variables (\(F, F', G, G'\)), and \(\mathrm{Re}\) as arguments. Boundary conditions are enforced in a dedicated function that reflects the physical constraints at the channel walls and centerline.

\subsection{Pressure Distribution}
A function computes the pressure coefficient \(c_p\) using \(\eta\), \(\mathrm{Re}\), and \(\alpha\) as in Equation~\ref{eq:pressure_coefficient}.

\subsection{Solving the BVP}
For each Reynolds number, an initial guess is provided to \texttt{solve\_bvp}. After convergence, velocity, pressure, stresses, and temperature profiles are computed and visualized.

\subsection{Visualization}
Five plots are produced for varying Reynolds numbers: velocity profile, pressure distribution, shear stress, normal stress in \(r\) and \(\theta\) directions, and the dimensionless temperature \(G\), providing a comprehensive view of the flow.

\section{Discussion and Conclusion}
This section discusses the results for velocity profiles, pressure, shear and normal stresses, and temperature distributions across several Reynolds numbers.

\subsection{Velocity Profile}
Figure~\ref{fig:velocity_profile} shows the dimensionless velocity \(F(\eta)\) versus \(\eta\). At \(\eta=0\) (wedge vertex), the velocity equals the free stream value (\(F'(0) = 1\)); as \(\eta \to 1\), the velocity drops to zero. Higher Reynolds numbers yield steeper gradients, indicating thinner boundary layers.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.8\textwidth]{figs/1.png}
\caption{Velocity profile for different Reynolds numbers.}
\label{fig:velocity_profile}
\vspace{10pt}
\end{figure}

The results agree with those reported by Panton~\cite{doi:https://doi.org/10.1002/9781118713075.ch14}, as shown in Figure~\ref{fig:velocity_comparison}, confirming the numerical method.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.8\textwidth]{figs/11.png}
\caption{Velocity ratio for positive and negative Reynolds numbers, compared with Panton~\cite{doi:https://doi.org/10.1002/9781118713075.ch14}.}
\label{fig:velocity_comparison}
\vspace{10pt}
\end{figure}

\subsection{Pressure Distribution}
Figure~\ref{fig:pressure_distribution} shows the pressure distribution along the centerline. For positive Reynolds numbers, pressure is positive; for negative, negative---indicating a favorable pressure gradient. Pressure increases with Reynolds number, reflecting greater momentum.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.8\textwidth]{figs/2.png}
\caption{Pressure distribution for different Reynolds numbers.}
\label{fig:pressure_distribution}
\vspace{10pt}
\end{figure}

At \(\eta=0\), pressure is near \(p_\infty + 0.61\); with increasing \(\eta\), it rises to a maximum before decreasing. Higher Reynolds numbers raise the pressure, while wider angles lower the gradient by reducing velocity gradients.

\subsection{Shear Stress}
Figure~\ref{fig:shear_stress} presents shear stress \(\tau_{\theta r}\) on the channel walls. Shear stress increases with Reynolds number, indicating stronger viscous effects at higher velocities.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.8\textwidth]{figs/3.png}
\caption{Shear stress for different Reynolds numbers.}
\label{fig:shear_stress}
\vspace{10pt}
\end{figure}

At \(\eta=0\), shear is zero (no velocity gradient at the wedge vertex). As \(\eta\) increases, shear rises, peaking for \(\mathrm{Re}=5\) due to steeper gradients near the wall.

\subsection{Normal Stress}
Figures~\ref{fig:normal_stress_r} and~\ref{fig:normal_stress_theta} show normal stresses \(\tau_{rr}\) and \(\tau_{\theta\theta}\) in the \(r\) and \(\theta\) directions, respectively. These are crucial for force and deformation predictions.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\begin{subfigure}[b]{0.45\textwidth}
    \centering
    \includegraphics[width=\textwidth]{figs/4.png}
    \caption{Normal stress in the \(r\)-direction.}
    \label{fig:normal_stress_r}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.45\textwidth}
    \centering
    \includegraphics[width=\textwidth]{figs/5.png}
    \caption{Normal stress in the \(\theta\)-direction.}
    \label{fig:normal_stress_theta}
\end{subfigure}
\caption{Normal stress components for different Reynolds numbers.}
\label{fig:normal_stress}
\vspace{10pt}
\end{figure}

Normal stresses increase with Reynolds number and decrease with larger angles. As \(\eta \to 1\), normal stresses approach zero, consistent with the boundary conditions.

\subsection{Dimensionless Temperature}
Figures~\ref{fig:temperature_with_dissipation} and~\ref{fig:temperature_no_dissipation} show the dimensionless temperature \(G\). Wall temperature is fixed at 300~K.

\subsubsection{With Viscous Dissipation}
Figure~\ref{fig:temperature_with_dissipation} displays \(G\) when viscous dissipation is included. At \(\eta=0\), \(G=0\); as \(\eta\) increases, \(G\) rises, then drops to zero at the boundary. Higher Reynolds numbers yield steeper temperature gradients, indicating thinner thermal boundary layers.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.8\textwidth]{figs/6.png}
\caption{Dimensionless temperature with viscous dissipation for different Reynolds numbers.}
\label{fig:temperature_with_dissipation}
\vspace{10pt}
\end{figure}

Larger channel angles produce flatter temperature profiles, i.e., thicker thermal boundary layers.

\subsubsection{Without Viscous Dissipation}
Figure~\ref{fig:temperature_no_dissipation} shows \(G\) with viscous heating neglected. In this case, the temperature is uniform (\(G=0\)) across the channel.

\begin{figure}[!htbp]
\centering
\vspace{10pt}
\includegraphics[width=0.8\textwidth]{figs/7.png}
\caption{Dimensionless temperature without viscous dissipation for different Reynolds numbers.}
\label{fig:temperature_no_dissipation}
\vspace{10pt}
\end{figure}

\subsection{Conclusion}
This numerical study analyzes Jeffery-Hamel flow using Python's \texttt{scipy.integrate.solve\_bvp}. The governing equations and boundary conditions were derived and solved for velocity, pressure, shear and normal stresses, and temperature distributions over a range of Reynolds numbers. Results are consistent with literature~\cite{doi:https://doi.org/10.1002/9781118713075.ch14}, confirming the approach. The final presentation is accessible at~\href{https://mega.nz/folder/tukxRRRZ#vAswA-ejytKRosnfE0HNgw}{this link}.

\newpage
{
\fontsize{12pt}{10pt}\selectfont
\bibliographystyle{ieeetr-fa}
\bibliography{etc/Panton}
\addcontentsline{toc}{section}{References}
}
\newpage
\end{document}
