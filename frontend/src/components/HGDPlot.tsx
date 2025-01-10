import { useEffect, useRef } from 'react';

interface HGDPlotProps {
  bins?: number;
  width?: number;
  height?: number;
}

export function HGDPlot({ bins = 30, width = 800, height = 600 }: HGDPlotProps) {
  const plotRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let ws: WebSocket | null = null;

    const connectWebSocket = () => {
      ws = new WebSocket('ws://localhost:65253/state?token=KlPNS0t9');

      ws.onopen = () => {
        console.log('Connected to httpgd WebSocket');
      };

      ws.onmessage = (event) => {
        const data = JSON.parse(event.data);
        console.log(data);

        if (data.type === 'render' && plotRef.current) {
          plotRef.current.innerHTML = data.content;
        }
      };

      ws.onerror = (error) => {
        console.error('WebSocket error:', error);
      };

      ws.onclose = () => {
        console.log('WebSocket closed, attempting to reconnect...');
        setTimeout(connectWebSocket, 1000);
      };
    };

    connectWebSocket();

    const updatePlot = async () => {
      try {
        await fetch('http://localhost:3838/plot', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ bins }),
        });
      } catch (error) {
        console.error('Error updating plot:', error);
      }
    };

    updatePlot();

    return () => {
      if (ws) {
        ws.close();
      }
    };
  }, [bins]);

  return (
    <div 
      ref={plotRef} 
      style={{ 
        width: `${width}px`, 
        height: `${height}px`,
        border: '1px solid #ccc',
        borderRadius: '4px',
        overflow: 'hidden'
      }} 
    />
  );
}
