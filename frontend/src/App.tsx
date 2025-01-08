import { useState, useEffect } from 'react';
import './App.css';

function App() {
  const [plotData, setPlotData] = useState<string | null>(null);
  const [message, setMessage] = useState<string>("");
  const [bins, setBins] = useState(30);

  useEffect(() => {
    // Set up event listeners
    const helloHandler = (event: CustomEvent) => {
      setMessage(event.detail.message);
    };

    const plotHandler = (event: CustomEvent) => {
      setPlotData(event.detail.plotData);
    };

    // Add event listeners
    window.addEventListener('helloResponse', helloHandler as EventListener);
    window.addEventListener('plotResponse', plotHandler as EventListener);

    // Fetch hello message
    fetch('http://localhost:3838/', {
      method: 'POST',
      credentials: 'include',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ getHello: true })
    });

    // Cleanup
    return () => {
      window.removeEventListener('helloResponse', helloHandler as EventListener);
      window.removeEventListener('plotResponse', plotHandler as EventListener);
    };
  }, []); // Empty dependency array means this runs once on mount

  // Fetch new plot when bins change
  useEffect(() => {
    fetch('http://localhost:3838/', {
      method: 'POST',
      credentials: 'include',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ getPlot: true, bins })
    });
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
        <div className="plot-container">
          <img 
            src={`data:image/png;base64,${plotData}`} 
            alt="R Plot" 
            style={{ maxWidth: '100%', height: 'auto' }}
          />
        </div>
      )}
    </div>
  );
}

export default App;
