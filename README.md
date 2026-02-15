# python_code_phd

Dieses Repository enthält beispielhaften Python-Code, den ich während meiner Promotion geschrieben habe (erstmals erstellt in 2022 und über die Jahre verfeinert, Funktionen hinzugefügt etc.). Der Code dient der Auswertung und schnellen Visualisierung von Messdaten aus der **Transienten Absorptionsspektroskopie (TA)**, also optical-pump-optical-probe. Der Datensatz einer TA-Messung beschreibt die Änderung der Absorption der Probe ΔA abhängig von der Abfragewellenlänge λ und der Verzögerungszeit t zwischen Anregungs- und Abfragepuls, also ΔA(λ,t). 

**ta.py** basiert auf dem Prinzip der Objekt-orientierten Programmierung (OOP). Die Klasse *TAData* wird definiert. Eine Instanz dieser Klasse beinhaltet eine TA-Messung, bei der Erzeugung einer Instanz werden also Daten sowie Metadaten einer Messung eingelesen und als Attribute der Instanz gespeichert. 

Mittels der Methoden von *TAData* können nun Korrekturen an den Daten vorgenommen werden, schnelle Visualisierungen in einheitlichem Design erzeugt werden (z.B. *plot_spectra* zur Darstellung differentieller Spektren zu definierten Verzögerungszeiten ΔA(λ)) oder aufwendigere Analysen (z.B. *analyse_bleach_shift* zur Untersuchung eines zeitlichen Shifts eines Signals, hier des sogenannten Bleach-Signals) durchgeführt werden.

Weiter unten wird die Klasse *TAMeasurementRow* definiert. Eine Instanz dieser Klasse beinhaltet eine Messreihe, also z.B. mehrere Einzelmessungen an einer Probe mit unterschiedlichen Anregungsbedingungen. Auch diese Klasse beinhaltet Methoden zur einfachen Visualisierung von Daten (z.B. *compare_spectra* zum Vergleich differentieller Anregungsspektren zu definierten Verzögerungszeiten ΔA(λ) unter unterschiedlichen Anregungsbedingungen) sowie Methoden, die aufwendigere Analysen (z.B. *mean_bleach_at_long_delays*) durchführen. 

**How_to_use_ta_py.ipynb** habe ich damals erstellt, um Studierenden, die den Code nutzen wollten, die Funktionsweise zu erläutern (inklusive absichtlich erzeugter Errors usw.) und dient der Veranschaulichung. Das File wurde allerdings nie vervollständigt, da ich mich dann in der Regel einfach einmal mit den Studis hingesetzt habe, um ihnen die Funktionsweise des Codes direkt zu zeigen.    


Im Ordner **further_exemplary_code** sind weitere Code-Beispiele hinterlegt, die ähnlich aufgebaut sind. Diese dienen der Auswertung von Messungen aus Terahertz-Experimenten (optp.py, **o**ptical-**p**ump-**t**erahertz-**p**robe) und von Raman-Messungen.
