\begin{tikzpicture}
    \begin{axis}[
        title={},
        width=0.50\textwidth,
        height=0.50\textwidth,
        xlabel={$\log(\text{slots})$},
        ylabel={$\log(1/\epsilon)$},
        xmin=3, xmax=16,
        ymin=8, ymax=40,
        xtick distance = 2,
        ytick distance = 2,
        legend pos=north west,
        ymajorgrids=true,
        grid style=dashed],
    ]
    \addplot+[color=red, mark=triangle*, error bars/.cd, y dir=both, y explicit]
        coordinates{
        {{.DataReal}}};
    \addplot+[color=blue, mark=square*, error bars/.cd, y dir=both, y explicit]
        coordinates{
        {{.DataImag}}
        };
    \legend{Real, Imag}
        
\end{axis}
\end{tikzpicture}