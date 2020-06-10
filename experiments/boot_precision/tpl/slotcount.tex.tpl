\begin{tikzpicture}
    \begin{axis}[
        title={},
        width=0.55\textwidth,
        height=0.55\textwidth,
        xlabel={$\log(\text{slots})$},
        ylabel={Precision (bits)},
        xmin=2, xmax=16,
        ymin=10, ymax=40,
        xtick=data,
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