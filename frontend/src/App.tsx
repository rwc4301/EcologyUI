import { useState, useEffect, createElement } from 'react';
import './App.css';

function App() {
  const [plotData, setPlotData] = useState<HTMLElement | null>(null);
  const [message, setMessage] = useState<string>("");
  const [bins, setBins] = useState(30);

  // Fetch hello message
  useEffect(() => {
    const fetchHello = async () => {
      try {
        const response = await fetch('http://localhost:3838/hello');
        console.log(response);
        const data = await response.json();
        setMessage(data.message);
      } catch (error) {
        console.error('Error fetching hello:', error);
        setMessage("Error connecting to R server");
      }
    };

    fetchHello();
  }, []);

  // Fetch plot when bins change
  useEffect(() => {
    const fetchPlot = async () => {
      try {
        const response = await fetch('http://localhost:3838/plot', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ bins })
        });
        const data = await response.json();
        setPlotData(data.plotData);
      } catch (error) {
        console.error('Error fetching plot:', error);
      }
    };

    fetchPlot();
  }, [bins]);

  return (
    <div className="App">
      <h1>Ecology Analyses</h1>
      
      <div className="message">
        {message || "Waiting for R server..."}
      </div>

      <div>
        <label>
          Number of bins:
          <input
            type="number"
            value={bins}
            onChange={(e) => setBins(Number(e.target.value))}
            min="1"
            max="100"
          />
        </label>
      </div>

      {plotData && (
        <div 
          className="plot-container"
          dangerouslySetInnerHTML={{ __html: plotData }}
        >
        </div>
      )}
    </div>
  );
}

export default App;
