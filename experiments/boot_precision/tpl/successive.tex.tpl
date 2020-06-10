\begin{tikzpicture}
    \begin{axis}[
        title={},
        width=1\textwidth,
        height=0.55\textwidth,
        xlabel={Iteration},
        ylabel={Precision (bits)},
        xmin=0, xmax=51,
        ymin=10, ymax=36,
        xtick distance = 5,
        ytick distance = 2,
        legend pos=north west,
        ymajorgrids=true,
        grid style=dashed],
    ]
    \addplot+[color=red, mark=triangle*, error bars/.cd, y dir=both, y explicit] coordinates {
        % Real
{{.DataReal}}
    };
    \addplot+[color=blue, mark=square*, error bars/.cd, y dir=both, y explicit] coordinates {
        % Imag
{{.DataImag}}
    };
    
    \legend{Real, Imag}
        
    \end{axis}
\end{tikzpicture}
