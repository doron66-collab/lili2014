import { useState } from 'react';
import NSCLCViewer from './NSCLCViewer';
import DataPanel from './DataPanel';
import IntroScreen from './IntroScreen';

function App() {
  const [started, setStarted] = useState(false);

  return (
    <div style={{ width: '100vw', height: '100vh', overflow: 'hidden', position: 'relative' }}>
      {!started && <IntroScreen onEnter={() => setStarted(true)} />}
      <NSCLCViewer />
      <DataPanel />
    </div>
  );
}

export default App;
